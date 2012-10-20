#!/usr/bin/env python

'''
\n\nThis script is INCOMPLETE and UNTESTED.\n\n
Grabs the contents of "snap_xaui" for analysis.
'''

#4/9/2011 - fixed antenna bug. Things seem to be working fine.

import corr, time, numpy, pylab, struct, sys, logging, pylibmc


#brams
brams=['bram_msb','bram_lsb','bram_oob']

# OOB signalling bit allocations:
linkdn_bit = 8
mrst_bit = 4
adc_bit = 3
eof_bit = 2
sync_bit = 1
hdr_bit = 0

def exit_fail():
    print 'FAILURE DETECTED. Log entries:\n',lh.printMessages()
    print "Unexpected error:", sys.exc_info()

    try:
        c.disconnect_all()
    except: pass
    raise
    exit()

def exit_clean():
    try:
        c.disconnect_all()
    except: pass
    exit()

def xaui_feng_unpack(xeng,xaui_port,bram_dump,hdr_index,pkt_len,skip_indices,mcache):
    pkt_64bit = struct.unpack('>Q',bram_dmp['bram_msb'][(4*hdr_index):(4*hdr_index)+4]+bram_dmp['bram_lsb'][(4*hdr_index):(4*hdr_index)+4])[0]
    pkt_mcnt =(pkt_64bit&((2**64)-(2**16)))>>16
    #pkt_ant  = xeng*c.config['n_xaui_ports_per_fpga']*c.config['n_ants_per_xaui'] + xaui_port*c.config['n_ants_per_xaui'] + pkt_64bit&((2**16)-1)
    #antenna variable in antenna design is already the actual antenna number. Dont neet to correct for it.
    pkt_ant  = pkt_64bit&((2**16)-1) 
    #print pkt_ant
    pkt_freq = pkt_mcnt%n_chans
    sum_polQ_r = 0
    sum_polQ_i = 0
    sum_polI_r = 0
    sum_polI_i = 0

    #average the packet contents - ignore first entry (header)
    for pkt_index in range(1,(pkt_len)):
        abs_index = hdr_index + pkt_index

        if skip_indices.count(abs_index)>0:
            #print 'Skipped %i'%abs_index
            continue

        pkt_64bit = struct.unpack('>Q',bram_dmp['bram_msb'][(4*abs_index):(4*abs_index)+4]+bram_dmp['bram_lsb'][(4*abs_index):(4*abs_index)+4])[0]
        raw_xaui_data = bram_dmp['bram_msb'][(4*abs_index):(4*abs_index)+4]+bram_dmp['bram_lsb'][(4*abs_index):(4*abs_index)+4]
        #save raw data to redis with the following format:
        #px?:snap_xaui_raw:antenna:channel = string of length 128 (each sample of a single channel ) * 1(change format to one byte instead of nibble) * 2(r,i) * 2(dual pol)
        if mcache.get('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq)) == None:
            mcache.set('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq), '')
        if len(mcache.get('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq))) == 256:
            #print 'writing Ant%d, Chan%d into redis.'%(pkt_ant,pkt_freq)
            mcache.set('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq), '')
            
   
        mcache.set('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq),mcache.get('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq)) + raw_xaui_data) 
        #if (pkt_ant==0 and pkt_freq == 650): 
            #print len(mcache.get('px%d:snap_xaui_raw:%d:%d'%(xeng+1,pkt_ant,pkt_freq)))
        for offset in range(48,-16,-16):
            polQ_r = (pkt_64bit >> (offset + 12)) & 0xf
            polQ_i = (pkt_64bit >> (offset +  8)) & 0xf
            polI_r = (pkt_64bit >> (offset +  4)) & 0xf
            polI_i = (pkt_64bit >> (offset     )) & 0xf
        
        
            #square each number and then sum it
            sum_polQ_r += (float(((numpy.int8(polQ_r << 4)>> 4)))/(2**binary_point))**2
            sum_polQ_i += (float(((numpy.int8(polQ_i << 4)>> 4)))/(2**binary_point))**2
            sum_polI_r += (float(((numpy.int8(polI_r << 4)>> 4)))/(2**binary_point))**2
            sum_polI_i += (float(((numpy.int8(polI_i << 4)>> 4)))/(2**binary_point))**2

        #print 'Processed %i. Sum Qr now %f...'%(abs_index, sum_polQ_r)

    num_accs = (pkt_len-len(skip_indices))*(64/16)

    level_polQ_r = numpy.sqrt(float(sum_polQ_r)/ num_accs)
    level_polQ_i = numpy.sqrt(float(sum_polQ_i)/ num_accs)
    level_polI_r = numpy.sqrt(float(sum_polI_r)/ num_accs)
    level_polI_i = numpy.sqrt(float(sum_polI_i)/ num_accs)

    rms_polQ = numpy.sqrt(((level_polQ_r)**2)  +  ((level_polQ_i)**2))
    rms_polI = numpy.sqrt(((level_polI_r)**2)  +  ((level_polI_i)**2))

    return {'pkt_mcnt': pkt_mcnt,\
            'pkt_ant':pkt_ant,\
            'pkt_freq':pkt_freq,\
            'rms_polQ':rms_polQ,\
            'rms_polI':rms_polI,\
            'level_polQ_r':level_polQ_r,\
            'level_polQ_i':level_polQ_i,\
            'level_polI_r':level_polI_r,\
            'level_polI_i':level_polI_i}

ant=0
if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('snap_xaui_spectrum_construct.py [options] CONFIG_FILE')
    p.set_description(__doc__)
    p.add_option('-a', '--ant', dest='ant', type='int',default=0,
        help="Select which antenna's data to plot. Default: 0")
    p.add_option('-v', '--verbose', dest='verbose', action='store_true',
        help='Be verbose about packet decoding.')
    p.add_option('-p', '--paranoid', dest='paranoid', action='store_true',
        help='Print paranoid raw output.')

    opts, args = p.parse_args(sys.argv[1:])

    if args==[]:
        print 'Please specify a configuration file! \nExiting.'
        exit()

    ant=opts.ant
    print ant

lh=corr.log_handlers.DebugLogHandler()
lh.setLevel(10)
try:
    print 'Parsing configuration file...',
    c=corr.corr_functions.Correlator(args[0],lh)
    for s,server in enumerate(c.config['servers']): c.loggers[s].setLevel(10)
    print 'done.'
    
    mcache = pylibmc.Client(['%s'%c.config['rx_udp_ip_str']])

    packet_len=c.config['10gbe_pkt_len']
    n_chans=c.config['n_chans']
    binary_point = c.config['feng_fix_pnt_pos']
    num_bits = c.config['feng_bits']
    adc_levels_acc_len = c.config['adc_levels_acc_len']

    if num_bits != 4:
        print 'This script is only written to work with 4 bit quantised values.'
        exit_clean()
    
    #figure out on which XAUI port this antenna resides:
    target_fpga, target_xaui, xaui_ant = c.get_ant_location(ant)

    dev_name = 'snap_xaui%i'%target_xaui

    print 'Looking for antenna %i on board %i (%s) XAUI port %i. Antenna %i on that port.'%(ant,target_fpga, c.servers[target_fpga],target_xaui,xaui_ant)

   # snap_rms_i= numpy.zeros(n_chans)
   # snap_rms_q= numpy.zeros(n_chans)
    snaps = [ numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans),
             numpy.zeros(n_chans)]
        
    raw_data_dump = 0

    offset=0 #channel offset to capture next
    pkt_len=0

    while offset < (n_chans-1):

        print 'Capturing frequency offset %i. Triggering at word offset %i...'%(offset,offset*pkt_len*c.config['n_ants_per_xaui']),
        sys.stdout.flush()
        bram_dmp=c.snap_x_single(c.fpgas[target_fpga],dev_name,brams,man_trig=False,wait_period=3,offset=offset*pkt_len*c.config['n_ants_per_xaui'])
        print 'got back data at offset %i.'%bram_dmp['offset']

        print 'Unpacking out of band contents...',
        sys.stdout.flush()
        bram_oob=dict()
        bram_oob={'raw':struct.unpack('>%iL'%(bram_dmp['length']),bram_dmp['bram_oob'])}
        #print bram_oob[f]['raw']
        bram_oob.update({'linkdn':[bool(i&(2**linkdn_bit)) for i in bram_oob['raw']]})
        bram_oob.update({'mrst':[bool(i&(2**mrst_bit)) for i in bram_oob['raw']]})
        bram_oob.update({'adc':[bool(i & (2**adc_bit)) for i in bram_oob['raw']]})
        bram_oob.update({'eof':[bool(i & (2**eof_bit)) for i in bram_oob['raw']]})
        bram_oob.update({'sync':[bool(i&(2**sync_bit)) for i in bram_oob['raw']]})
        bram_oob.update({'hdr':[bool(i & (2**hdr_bit)) for i in bram_oob['raw']]})
        print 'Done.'

        print 'Unpacking data contents...'
        skip_indices=[]
        pkt_hdr_idx = -1
        for i in range(0,bram_dmp['length']):
            if opts.paranoid:
                pkt_64bit = struct.unpack('>Q',bram_dmp['bram_msb'][(4*i):(4*i)+4]+bram_dmp['bram_lsb'][(4*i):(4*i)+4])[0]
                print '[%4i]: %016X'%(i,pkt_64bit),
                if bram_oob['eof'][i]: print '[EOF]',
                if bram_oob['linkdn'][i]: print '[LINK DOWN]',
                if bram_oob['mrst'][i]: print '[MRST]',
                if bram_oob['adc'][i]: print '[ADC_UPDATE]',
                if bram_oob['sync'][i]: print '[SYNC]',
                if bram_oob['hdr'][i]: print '[HDR]',
                print ''


            if bram_oob['linkdn'][i]:
                print 'LINK DOWN AT %i'%(i)
            elif bram_oob['adc'][i]:
                adc,amp = struct.unpack('>II',bram_dmp['bram_msb'][(4*i):(4*i)+4]+bram_dmp['bram_lsb'][(4*i):(4*i)+4])
                adc_ant,adc_pol=c.get_ant_index(target_fpga,target_xaui,adc&0xff)
                print ' ADC amplitude update at index %i for input %i (ant %i,%s) with RMS count of %3f.'%(i,adc,adc_ant,adc_pol,numpy.sqrt(float(amp)/adc_levels_acc_len))

                skip_indices.append(i) #skip_indices records positions in table which are ADC updates and should not be counted towards standard data.
            elif bram_oob['hdr'][i]:
                pkt_hdr_idx = i
                skip_indices=[]
    #            print ('HEADER RECEIVED')
            elif bram_oob['eof'][i]:
                #skip the first packet entry which has no header (snap block triggered on sync)
                if pkt_hdr_idx<0: continue
    #            print ('EOF RECEIVED')
                #packet_len is length of data, not including header, in 64-bit words
                pkt_len=i-pkt_hdr_idx+1
                if pkt_len-len(skip_indices) != (packet_len+1):
                    print 'MALFORMED PACKET! of length %i starting at index %i'%(pkt_len-len(skip_indices),i)
                    print 'skip_indices: (%i numbers):'%len(skip_indices),skip_indices
                else:
                    feng_unpkd_pkt=xaui_feng_unpack(target_fpga,target_xaui,bram_dmp,pkt_hdr_idx,pkt_len,skip_indices,mcache)
                    if opts.verbose: print '[Pkt@ %4i Len: %2i]     (MCNT %12u ANT: %1i, Freq: %4i)    {4 bit: Qr: %1.2f Qi: %1.2f Ir %1.2f Ii: %1.2f}'%(\
                        pkt_hdr_idx,\
                        pkt_len-len(skip_indices),\
                        feng_unpkd_pkt['pkt_mcnt'],\
                        feng_unpkd_pkt['pkt_ant'],\
                        feng_unpkd_pkt['pkt_mcnt']%n_chans,\
                        feng_unpkd_pkt['level_polQ_r'],\
                        feng_unpkd_pkt['level_polQ_i'],\
                        feng_unpkd_pkt['level_polI_r'],\
                        feng_unpkd_pkt['level_polI_i'])


                    if feng_unpkd_pkt['pkt_ant']%4 == 0:  
                        snaps[0][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polQ']
                        snaps[1][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polI']
                    if feng_unpkd_pkt['pkt_ant']%4 == 1:  
                        snaps[2][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polQ']
                        snaps[3][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polI']
                    if feng_unpkd_pkt['pkt_ant']%4 == 2:  
                        snaps[4][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polQ']
                        snaps[5][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polI']
                    if feng_unpkd_pkt['pkt_ant']%4 == 3:  
                        snaps[6][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polQ']
                        snaps[7][feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polI']

                    #if feng_unpkd_pkt['pkt_ant'] == ant:
                    #    snap_rms_q[feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polQ']
                    #    snap_rms_i[feng_unpkd_pkt['pkt_freq']]=feng_unpkd_pkt['rms_polI']

                    #if (offset>(n_chans-1)): 
                    #    #Got all the channels, exit.
                    #    offset=feng_unpkd_pkt['pkt_freq']
                    #    break
                    if (offset<n_chans):
                        if feng_unpkd_pkt['pkt_ant']%4 == 3:
                            offset=feng_unpkd_pkt['pkt_freq']+1
                            if offset == 2048: 
                                print 'got last channel'
                                break
                        else: offset = feng_unpkd_pkt['pkt_freq']
                    elif (offset-2 > feng_unpkd_pkt['pkt_freq']):
                        print 'Snap block failed to capture at the correct trigger offset. Exiting.'
                        exit_clean()
                    print offset
    print '\n\nDone Capturing. Plotting...' 
    
    for input in range(c.config['n_ants_per_xaui']):
        mcache.set('snap_xaui:%d'%((ant/4)*4*2+input),snaps[input].tostring())

    ant2plt = 2*(ant%4)
    pylab.figure(ant)
    ax1=pylab.subplot(211)
    pylab.title('Antenna %i\n"X" input rms. RMS = %f'%(ant,numpy.average(snaps[ant2plt])))
    pylab.plot(range(0,n_chans),snaps[ant2plt])
    #pylab.plot(range(0,n_chans),snap_rms_q,'.')
    pylab.setp(ax1.get_xticklabels(), visible=False)

    pylab.subplot(212,sharex=ax1,sharey=ax1)
    pylab.title('"Y" input rms. RMS = %f'%(numpy.average(snaps[ant2plt+1])))
    #pylab.plot(range(0,n_chans),snap_rms_i,'.',label='ant%i'%xaui_ant)
    pylab.plot(range(0,n_chans),snaps[ant2plt +1])
    pylab.ylim((0,1))
    pylab.xlim(0,n_chans)
    pylab.xlabel('Channel')
    
   # pylab.figure(ant+1)
   # ax1=pylab.subplot(211)
   # pylab.title('Antenna %i\n"X" input rms'%ant)
   # pylab.plot(range(0,n_chans),snap_rms_q[::-1])
   # #pylab.plot(range(0,n_chans),snap_rms_q,'.')
   # pylab.setp(ax1.get_xticklabels(), visible=False)

   # pylab.subplot(212,sharex=ax1,sharey=ax1)
   # pylab.title('"Y" input rms')
   # #pylab.plot(range(0,n_chans),snap_rms_i,'.',label='ant%i'%xaui_ant)
   # pylab.plot(range(0,n_chans),snap_rms_i[::-1])
   # pylab.ylim((0,1))
   # pylab.xlim(0,n_chans)
   # pylab.xlabel('Channel')

    pylab.show()

except SystemExit:
    exit_clean()

except KeyboardInterrupt:
    lh.printMessages()
    exit_clean()
except:
    exit_fail()

exit_clean()


