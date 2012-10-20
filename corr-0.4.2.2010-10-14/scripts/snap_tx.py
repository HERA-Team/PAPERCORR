#! /usr/bin/env python
import corr, os, numpy, sys, optparse, time,struct

regs = ['bram_msb', 'bram_lsb', 'bram_oob']
basename = 'switch_txsnap'

#OOB signaling bit offsets:
ip_addr_bit_width = 32-8
ip_addr_bit_offset = 5
eof_bit = 4
link_up_bit = 3
tx_led_bit = 2
tx_afull_bit = 1
tx_over_bit = 0

#NSAMP = 1024
NSAMP = 2048

def ip2str(pkt_ip):
    ip_4 = (pkt_ip&((2**32)-(2**24)))>>24
    ip_3 = (pkt_ip&((2**24)-(2**16)))>>16
    ip_2 = (pkt_ip&((2**16)-(2**8)))>>8
    ip_1 = (pkt_ip&((2**8)-(2**0)))>>0
    #print 'IP:%i. decoded to: %i.%i.%i.%i'%(pkt_ip,ip_4,ip_3,ip_2,ip_1)
    return '%i.%i.%i.%i'%(ip_4,ip_3,ip_2,ip_1)    


def snap_x_single(fpga,dev_name,brams,man_trig=False,man_valid=False,wait_period=1,offset=-1):
    """Grabs from a single X engine board in the system, brams from device.\n
    \tfpga: katcp_wrapper.fpga object from which to grab snapshot data.\n
    \tdev_name: string, name of the snap block.\n
    \tman_trig: boolean, Trigger the snap block manually.\n
    \toffset: integer, wait this number of valids before beginning capture. Set to negative value if your hardware doesn't support this or the circular capture function.\n
    \tcircular_capture: boolean, Enable the circular capture function.\n
    \twait_period: integer, wait this number of seconds between triggering and trying to read-back the data.\n
    \tbrams: list, names of the bram components.\n
    \tRETURNS: dictionary with keywords: \n
    \t\tlengths: list of integers matching number of valids captured off each fpga.\n
    \t\toffset: optional (depending on snap block version) list of number of valids elapsed since last trigger on each fpga.
    \t\t{brams}: list of data from each fpga for corresponding bram.\n"""
    #2010-02-19 JRM Updated to match snap_x.
    #WORKING OK 2009-07-01

    #print 'Triggering Capture...',
    fpga.write_int(dev_name+'_ctrl',(0 + (man_trig<<1) + (man_valid<<2)))
    fpga.write_int(dev_name+'_ctrl',(1 + (man_trig<<1) + (man_valid<<2)))

    time.sleep(wait_period)

    done=False
    start_time=time.time()
    while (not done) or ((time.time() - start_time) < wait_period):
        addr = fpga.read_uint(dev_name+'_addr')
        if addr == NSAMP-1 : done = True 
        #if addr == 1023 : done = True 
        else : done = False
    bram_dmp=dict()

    bram_size= fpga.read_uint(dev_name+'_addr')&0x7fffffff
    bram_dmp={'length':bram_size+1}
    bram_dmp['offset']=0
    if (bram_size != fpga.read_uint(dev_name+'_addr')&0x7fffffff) or bram_size==0:
        #if address is still changing, then the snap block didn't finish capturing. we return empty.
        print "Looks like snap block didn't finish."
        bram_dmp['length']=0
        bram_dmp['offset']=0
        bram_size=0

    if (bram_dmp['offset'] < 0):
        raise RuntimeError('SNAP block hardware or logic failure happened. Returning no data.')
    #if offset<0, there was a big error. this should never happen. unless you held stop high while resetting before inputting valid data during a circular capture? In any case, zero the output.
        bram_dmp['length']=0
        bram_dmp['offset']=0
        bram_size=0

    for b,bram in enumerate(brams):
        bram_path = dev_name+'_'+bram
        if (bram_size == 0):
            bram_dmp[bram]=[]
        else:
            bram_dmp[bram]=(fpga.read(bram_path,(bram_size+1)*4))
    return bram_dmp


try:

    import optparse 
    o = optparse.OptionParser()
    o.add_option('-f','--fromw', dest='fromw', type='string', default='switch',
            help='Where to snap this data from the switch tx snap block pr gpu tx snap block. "gpu" has 2048 data points per fpga where as witch only has 1024.keywords are "gpu" or "switch"')
    o.add_option('-r', '--roach', dest='roach', type=int, 
            help='specify a specific roach to capture/report from.')
    opts, args = o.parse_args(sys.argv[1:])


    if args==[]:
        print 'Please specify a configuration file! \n Exiting.'
        exit()

    lh = corr.log_handlers.DebugLogHandler()
    p = corr.corr_functions.Correlator(args[0], lh)
    for s, server in enumerate(p.config['servers']): p.loggers[s].setLevel(10)

    print opts.roach

    if opts.roach != None:
        servers = [p.fpgas[opts.roach]]
    else : servers = p.fpgas

    if opts.fromw == 'switch':
        basename = 'switch_txsnap'
        NSAMP=1024
    elif opts.fromw == 'gpu':
        basename='gpu_txsnap'
        NSAMP=2048
    

    report = dict()
    binary_point = 3
    num_bits = 4 
    packet_len = 1024
    n_ants_per_roach = 4 
    
    print '------------------------'
    print 'Grabbing snap data...'
    for f, fpga in enumerate(servers):
        report['px%i:bram_dmp'%f] = snap_x_single(fpga, dev_name = basename, brams = regs, man_trig=False, man_valid=False, wait_period=3)
    print 'done'

    print 'Unpacking bram contents....',
    sys.stdout.flush()
    bram_oob = dict()
    bram_data = dict()

    for f, fpga in enumerate(servers):
        bram_oob['px%d'%f] = {'raw':struct.unpack('>%iL'%NSAMP, report['px%d:bram_dmp'%f]['bram_oob'])}
        bram_data['px%d'%f] = {'rawmsb':struct.unpack('>%iL'%NSAMP, report['px%d:bram_dmp'%f]['bram_msb'])}
        bram_data['px%d'%f].update({'rawlsb':struct.unpack('>%iL'%NSAMP, report['px%d:bram_dmp'%f]['bram_lsb'])})
        bram_oob['px%d'%f].update({'eof':[bool(i & (1<<eof_bit)) for i in bram_oob['px%d'%f]['raw']]})
        bram_oob['px%d'%f].update({'link':[bool(i & (1<<link_up_bit)) for i in bram_oob['px%d'%f]['raw']]})
        bram_oob['px%d'%f].update({'tx_led':[bool(i & (1<<tx_led_bit)) for i in bram_oob['px%d'%f]['raw']]})
        bram_oob['px%d'%f].update({'tx_afull':[bool(i & (1<<tx_afull_bit)) for i in bram_oob['px%d'%f]['raw']]})
        bram_oob['px%d'%f].update({'tx_over':[bool(i & (1<<tx_over_bit)) for i in bram_oob['px%d'%f]['raw']]})
        bram_oob['px%d'%f].update({'ip_addr':[(i>>5)  for i in bram_oob['px%d'%f]['raw']]})
        bram_data['px%d'%f].update({'dataword':[(bram_data['px%d'%f]['rawmsb'][i]<<32) + bram_data['px%d'%f]['rawlsb'][i] for i in range(NSAMP)]})
    print 'done unpacking.'

    print 'Analysing packets:'

    #for f,fpga in enumerate(servers):
    for f,fgpa in enumerate(servers):
        for addr in range(NSAMP):
            if opts.roach != None:
                print '[px%i] IDX: %i Contents:' %(opts.roach, addr),
            else:    
                print '[px%i] IDX: %i Contents:' %(f, addr),
            print '[%s]'%ip2str(bram_oob['px%d'%f]['ip_addr'][addr]),
            if bram_oob['px%d'%f]['link'][addr]: print '[link] ',
            if bram_oob['px%d'%f]['eof'][addr]: print '[eof] ',
            if bram_oob['px%d'%f]['tx_led'][addr]: print '[tx_led] ',
            if bram_oob['px%d'%f]['tx_afull'][addr]: print '[tx_afull] ',
            if bram_oob['px%d'%f]['tx_over'][addr]: print '[tx_over] ',
            print 'DATUM = ', hex(bram_data['px%d'%f]['dataword'][addr]),
            print ''


         


except(KeyboardInterrupt):
    print 'Keyboard Interrupt'
    exit()


    

