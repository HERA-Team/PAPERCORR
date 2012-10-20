#!/usr/bin/env python

"""
Transmits correlator vacc data by reading out a snap block on the roach.

Author: Jason Manley
Date: 2010/01/20

Revisions:
2011-05-05  ZA  Added read out of corr_read_missing registers and saved 
                to memcachd.
    
2011-04-25  ZA  Start changes to use memcached python library to read out 
                corr_read missing stuff and adc amplitudes etc...

2010-02-24  JRM Added support for multiple X engines per FPGA.
                Now prints time difference between integration dumps.

2010-01-20  JRM Mangled together from other bits of code.
old stuff:
2008-09-10  JRM Bugfix "pack"

2008-02-13  JRM Further cleanups
                Additional sanity checks

2008-02-08  JRM New packet format
                Removed Rawpacket class - unified with CasperN_Packet
                Neatened CasperN_RX_correlator

2007-08-29  JRM Changed some endian-ness handling for packet decoding
"""

import time, os, socket, struct, sys, pylibmc

#  ____                          _   _     ____            _        _   
# / ___|__ _ ___ _ __   ___ _ __| \ | |   |  _ \ __ _  ___| | _____| |_ 
#| |   / _` / __| '_ \ / _ \ '__|  \| |   | |_) / _` |/ __| |/ / _ \ __|
#| |__| (_| \__ \ |_) |  __/ |  | |\  |   |  __/ (_| | (__|   <  __/ |_ 
# \____\__,_|___/ .__/ \___|_|  |_| \_|___|_|   \__,_|\___|_|\_\___|\__|
#               |_|                  |_____|                            


class CasperN_Packet:
    """Pack and unpack the binary correlation data in a CasperN packet,
    assuming the data is stored as signed (4 byte) integers."""

    def __init__(self, endian='>'):
        self.data_fmt = 'i'
        self.word_size = struct.calcsize(self.data_fmt)
        self.endian = endian
        self.header_fmt = '%sHHHHQQQQQQ' % (endian)
        self.header_size = struct.calcsize(self.header_fmt)


    def get_hdr_size(self):
        return self.header_size

    def pack_from_prms(self, timestamp, xeng, offset, heap_len, data):
        """Create a packet."""
        if type(data) is str: 
            str_data = data
        else:
            fmt = '%s%d%s' % (self.endian, len(d['data'])*self.word_size,self.data_fmt)
            str_data = struct.pack(fmt, 'data')

        option1 = (50<<48) + (3<<32) + (0<<16) + xeng
        option2 = (51<<48) + (len(str_data)<<24)
        option3 = (3<<48) + timestamp
        option4 = (4<<48) + heap_len
        option5 = (5<<48) + offset
        option6 = (1<<63) + (53<<48) + 0

        #print "PKT sending at timestamp %i for xeng %i at offset %i."%(timestamp,xeng,offset)

        return struct.pack(self.header_fmt,0x4b52,3,0,6,option1, option2, option3, option4, option5, option6) + str_data


#  _______  __    ____             _        _
# \_   _\ \/ /   / ___|  ___   ___| | _____| |_
#   | |  \  /    \___ \ / _ \ / __| |/ / _ \ __|
#   | |  /  \     ___) | (_) | (__|   <  __/ |_
#   |_| /_/\_\___|____/ \___/ \___|_|\_\___|\__|
#           |_____|

class UDP_TX_Socket(socket.socket):
    """Implements a UDP socket which transmits at a given ip, port."""
    def __init__(self, ip, port):
        self.ip = ip
        #print 'Set ip to %s' %self.ip
        self.port = port
        socket.socket.__init__(self, type=socket.SOCK_DGRAM)
        self.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        #self.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, 20)
        #self.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_LOOP, 1)
        #req = struct.pack('4sl', socket.inet_aton(ip), socket.INADDR_ANY)
        #self.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, req)
    def tx(self, data):
        """Send a UDP packet containing binary 'data'."""
        #print 'Sending to ip %s on port %i' %(self.ip,self.port)
        return self.sendto(data, (self.ip, self.port))

class CasperN_TX_Socket(UDP_TX_Socket, CasperN_Packet):
    """Combines a UDP_TX_Socket with the casper_n packet format."""
    def __init__(self, ip, port, endian='>'):
        UDP_TX_Socket.__init__(self, ip, port)
        CasperN_Packet.__init__(self, endian=endian)
    def send_packet(self, timestamp, xeng, offset, heap_len, data):
        """Send a UDP packet using the casper_n packet format."""
        return self.tx(self.pack_from_prms(timestamp, xeng, offset, heap_len, data))


class CorrTX:
    def __init__(self, pid, endian='>',ip='10.0.0.1', x_per_fpga=2, port=7147, payload_len=4096, verbose=False, timestamp_rnd=1024*128):
        self.pid=pid
        self.endian = endian
        self.casper_sock=CasperN_TX_Socket(ip,port,endian)
        self.ip=ip
        self.port=port
        self.payload_len=payload_len
        self.verbose=verbose
        self.x_per_fpga = x_per_fpga
        self.mcache = pylibmc.Client([ip])

        self.corr_read_missing = self.corr_read_missing_init(pid)

        self.ant_levels = open('/proc/%i/hw/ioreg/ant_levels'%(pid),'r')
        self.ant_levels_mean = open('/proc/%i/hw/ioreg/ant_levels_mean'%(pid),'r')
        self.adc_trigger = open('/proc/%i/hw/ioreg/adc_level_start'%(pid),'w')

        self.snap_xaui_trig_offset = open('/proc/%i/hw/ioreg/snap_xaui0_trig_offset'%(pid),'w')
        self.snap_xaui_ctrl = open('/proc/%i/hw/ioreg/snap_xaui0_ctrl'%(pid),'w')
        self.snap_xaui_oob = open('/proc/%i/hw/ioreg/snap_xaui0_bram_oob'%(pid), 'r')
        self.snap_xaui_msb = open('/proc/%i/hw/ioreg/snap_xaui0_bram_msb'%(pid), 'r')
        self.snap_xaui_lsb = open('/proc/%i/hw/ioreg/snap_xaui0_bram_lsb'%(pid), 'r')
        self.snap_xaui_addr = open('/proc/%i/hw/ioreg/snap_xaui0_addr'%(pid), 'r')
        self.pkt_len = 0
        self.freq_ant = {}
        for i in range(2048):
            self.freq_ant[i] = []
        self.finished_freq = []
        
        self.snap_addr=[]
        self.snap_bram=[]
        self.snap_en=[]
        self.xeng=[]
        self.vacc_mcnt_l=[]
        self.vacc_mcnt_h=[]
        for x in range(x_per_fpga):
            self.snap_addr.append(open('/proc/%i/hw/ioreg/snap_vacc%i_addr'%(pid,x),'r'))
            self.snap_bram.append(open('/proc/%i/hw/ioreg/snap_vacc%i_bram'%(pid,x),'r'))
            self.snap_en.append(open('/proc/%i/hw/ioreg/snap_vacc%i_ctrl'%(pid,x),'w'))
            self.vacc_mcnt_l.append(open('/proc/%i/hw/ioreg/vacc_mcnt_l%i'%(pid,x),'r'))
            self.vacc_mcnt_h.append(open('/proc/%i/hw/ioreg/vacc_mcnt_h%i'%(pid,x),'r'))
            #self.vacc_mcnt=(open('/proc/%i/hw/ioreg/vacc_mcnt%i'%(pid,x),'r'))

            xeng_file=(open('/proc/%i/hw/ioreg/inst_xeng_id%i'%(pid,x),'r'))
            xeng_file.seek(2)
            self.xeng.append(struct.unpack('>H',xeng_file.read(2))[0])
            xeng_file.close()
            print ('Ready to send output data from Xeng %i to IP %s on port %i.' %(self.xeng[x],ip,port))
    

        self.timestamp_rnd=timestamp_rnd
        self._tx()

    def corr_read_missing_init(self, pid, n_xaui = 1, x_per_fpga = 2):
        corr_read_missing = { 'xaui_errors' : [open('/proc/%i/hw/ioreg/xaui_err%i'%(pid,x),'r') for x in range(n_xaui)],
        'xaui_rx_cnt' : [open('/proc/%i/hw/ioreg/xaui_cnt%i'%(pid,x),'r') for x in range(n_xaui)],
        'gbe_tx_cnt' : [open('/proc/%i/hw/ioreg/gbe_tx_cnt%i'%(pid,x),'r') for x in range(n_xaui)],
        'gbe_tx_err' : [open('/proc/%i/hw/ioreg/gbe_tx_err_cnt%i'%(pid,x),'r') for x in range(n_xaui)],
        'rx_cnt' : [open('/proc/%i/hw/ioreg/rx_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'gbe_rx_cnt' : [open('/proc/%i/hw/ioreg/gbe_rx_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'gbe_rx_err_cnt' : [open('/proc/%i/hw/ioreg/gbe_rx_err_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'gbe_rx_down' : [open('/proc/%i/hw/ioreg/gbe_rx_down'%(pid),'r') for x in range(min(n_xaui, x_per_fpga))],
        'rx_err_cnt' : [open('/proc/%i/hw/ioreg/rx_err_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'loop_cnt' : [open('/proc/%i/hw/ioreg/loop_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'loop_error_cnt' : [open('/proc/%i/hw/ioreg/loop_err_cnt%i'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'mcnts' : [open('/proc/%i/hw/ioreg/loopback_mux%i_mcnt'%(pid,x),'r') for x in range(min(n_xaui,x_per_fpga))],
        'x_cnt' : [open('/proc/%i/hw/ioreg/pkt_reord_cnt%i'%(pid,x),'r') for x in range(x_per_fpga)],
        'x_miss' : [open('/proc/%i/hw/ioreg/pkt_reord_err%i'%(pid,x),'r') for x in range(x_per_fpga)],
        'last_miss_ant' : [open('/proc/%i/hw/ioreg/last_missing_ant%i'%(pid,x),'r') for x in range(x_per_fpga)],
        'vacc_cnt' : [open('/proc/%i/hw/ioreg/vacc_cnt%i'%(pid,x),'r') for x in range(x_per_fpga)],
        'vacc_err_cnt' : [open('/proc/%i/hw/ioreg/vacc_err_cnt%i'%(pid,x),'r') for x in range(x_per_fpga)]}
        return corr_read_missing

    def set_multi_ints_no_pickle(self, dict):
        #print dict
        for key in dict.keys():
            try:
                self.mcache.set(key, struct.pack(">" + "4s"*len(dict[key]), *dict[key]))
            except Exception, e: print 'REDIS ERROR (set multi ints no pickel)',e

    def get_corr_read_missing(self,corr_read_dictionary = {}): 
        corr_read2write = {}
        #print 'Reading corr_read_missing registers'
        for key in corr_read_dictionary.keys():
            corr_read2write['px%d:'%(self.xeng[0]+1)+key] = []
            for i in range(len(corr_read_dictionary[key])):
                corr_read_dictionary[key][i].flush()
                corr_read_dictionary[key][i].seek(0)
                corr_read_dictionary[key][i].flush()
                corr_read2write['px%d:'%(self.xeng[0]+1)+key].append(corr_read_dictionary[key][i].read())
        #print 'Saving corr_read_missing registers to redis'
        self.set_multi_ints_no_pickle(corr_read2write)     
        #print 'done'
        
    def adc_amplitudes(self):
        #4bytes for 1 input. There are 8 inputs.
        mem_size = 2 * 4 * 4   
        #print 'getting adc data,triggering adc level snap'
        self.adc_trigger.flush()
        self.adc_trigger.seek(0)
        self.adc_trigger.write(struct.pack('I',1))
        self.adc_trigger.flush()
        self.adc_trigger.seek(0)
        self.adc_trigger.write(struct.pack('I',0))
        self.adc_trigger.flush()
        time.sleep(.0005)
        self.ant_levels.flush()
        self.ant_levels.seek(0)
        self.ant_levels.flush()
        adc = self.ant_levels.read(mem_size)
        self.ant_levels_mean.flush()
        self.ant_levels_mean.seek(0)
        self.ant_levels_mean.flush()
        adc_mean = self.ant_levels_mean.read(mem_size)
        #print 'saving adc data into redis'
        try:
            self.mcache.set_multi({'px%d:adc_sum_squares'%(self.xeng[0]+1):adc, 'px%d:adc_sum'%(self.xeng[0]+1):adc_mean})
            #print 'px%d:adc_sum_squares'%(self.xeng[0]+1)
            #print 'done'
        except Exception, e: print 'REDIS ERROR (adc_amplitudes)',e

    def snap_xaui_ram(self,pkt_len,offset=-1, wait = 1):
        if offset >=0:
            #print 'freq offset in snap_xaui_ram=',offset,offset*32*4
            self.snap_xaui_trig_offset.flush()
            self.snap_xaui_trig_offset.seek(0)
            self.snap_xaui_trig_offset.write(struct.pack('I',offset*pkt_len*4))
            self.snap_xaui_trig_offset.flush()
            
        self.snap_xaui_ctrl.flush()
        self.snap_xaui_ctrl.seek(0)
        self.snap_xaui_ctrl.write(struct.pack('I',0))
        self.snap_xaui_ctrl.flush()
        self.snap_xaui_ctrl.seek(0)
        self.snap_xaui_ctrl.write(struct.pack('I',1))
        self.snap_xaui_ctrl.flush()
        time.sleep(wait)
        done = False
        start_time = time.time()
        while not (done and (offset > 0)) and ((time.time() - start_time) < wait):
            #print 'not done'
            self.snap_xaui_addr.seek(0)
            addr = struct.unpack('I',self.snap_xaui_addr.read())[0]
            #print addr
            done = bool(addr & 0x80000000)
        bram_dmp = dict()
        bram_size = (addr&0x7fffffff)
        bram_dmp = {'length':bram_size + 1}
        bram_dmp['offset'] = offset
        self.snap_xaui_addr.seek(0)
        if addr == struct.unpack('I',self.snap_xaui_addr.read())[0]:
            #print 'In data read loop'
            #begin read out of data.
            self.snap_xaui_oob.seek(0)
            bram_dmp['oob_data'] = self.snap_xaui_oob.read((bram_size+1)*4)
            self.snap_xaui_lsb.seek(0)
            bram_dmp['lsb_data'] = self.snap_xaui_lsb.read((bram_size+1)*4)
            self.snap_xaui_msb.seek(0)
            bram_dmp['msb_data'] = self.snap_xaui_msb.read((bram_size+1)*4)
        else: 
            self.snap_xaui_addr.seek(0)
            #print addr,struct.unpack('I',self.snap_xaui_addr.read())
        #print bram_dmp.keys()
        return bram_dmp 
        
    def xaui_unpack(self,bram_dmp, hdr_index,pkt_len,skip_indices,mcache):
        pkt_64bit_hdr = struct.unpack('Q', bram_dmp['msb_data'][(4*hdr_index):(4*hdr_index)+4] + bram_dmp['lsb_data'][(4*hdr_index):(4*hdr_index+4)])[0] 
        pkt_mcnt = pkt_64bit_hdr >>16
        pkt_ant = pkt_64bit_hdr & 0xffff
        pkt_freq = pkt_mcnt%2048 #mcnt % nchans
        raw_xaui_data=''
        for pkt_index in range(1,(pkt_len)):
            abs_index = hdr_index + pkt_index
            if skip_indices.count(abs_index)>0:continue

            raw_xaui_data += bram_dmp['msb_data'][(4*abs_index):(4*abs_index)+4]+bram_dmp['lsb_data'][(4*abs_index):(4*abs_index)+4]
            if len(raw_xaui_data) == 256:
                #print 'writing Ant%d, Chan%d into redis.'%(pkt_ant,pkt_freq)
                try:
                    mcache.set('px%d:snap_xaui_raw:%d:%d'%(self.xeng[0]+1,pkt_ant%4,pkt_freq), raw_xaui_data)
                except Exception, e: print 'REDIS ERROR(xaui unpack)',e
        return pkt_ant,pkt_freq



    def xaui_parse(self,bram_dmp):
        bram_oob = {'raw':struct.unpack('%iL'%(bram_dmp['length']), bram_dmp['oob_data'])}
        bram_oob.update({'linkdn':[bool(i&(1<<8)) for i in bram_oob['raw']]})
        bram_oob.update({'mrst':[bool(i&(1<<4)) for i in bram_oob['raw']]})
        bram_oob.update({'adc':[bool(i&(1<<3)) for i in bram_oob['raw']]})
        bram_oob.update({'eof':[bool(i&(1<<2)) for i in bram_oob['raw']]})
        bram_oob.update({'sync':[bool(i&(1<<1)) for i in bram_oob['raw']]})
        bram_oob.update({'hdr':[bool(i & (1<<0)) for i in bram_oob['raw']]})
        if opts.verbose:
            for k in range(len(bram_oob['raw'])):
                if bram_oob['linkdn'][k]:print k,'linkdn'
                if bram_oob['mrst'][k]:print k,'mrst'
                if bram_oob['adc'][k]:print k,'adc'
                if bram_oob['eof'][k]:print k,'eof'
                if bram_oob['sync'][k]:print k,'sync'
                if bram_oob['hdr'][k]:print k,'hdr'
        
        skip_indices = []
        pkt_hdr_idx = -1
        
        #print 'bram_dmp["length"] = ', bram_dmp['length']
        if bram_dmp["length"] == 1:
            return max(self.finished_freq),self.pkt_len
        else: 
            for i in range(bram_dmp['length']):
                if bram_oob['adc'][i]:
                    #print 'Adding skip indices',i
                    skip_indices.append(i)
                elif bram_oob['hdr'][i]:
                    #print 'header at',i
                    pkt_hdr_idx = i
                    skip_indices = []
                elif bram_oob['eof'][i]:
                    #print 'got eof'
                    if pkt_hdr_idx<0:continue
                    self.pkt_len = i-pkt_hdr_idx+1
                    #if pkt_len-len(skip_indices) != 33:
                    #    pass
                    #print 'unpacking data'
                    ant,freq = self.xaui_unpack(bram_dmp,pkt_hdr_idx,self.pkt_len,skip_indices,self.mcache)
                    self.freq_ant[freq].append(ant) 
            for k in self.freq_ant.keys():
                if len(self.freq_ant[k]) == 4:
                    if (k in self.finished_freq)==False:
                        self.finished_freq.append(k)
                        #print 'appending'
                   #     print self.finished_freq
                if max(self.finished_freq) == 2047:
                    print 'got all channels.Initializing finished_freq...'
                    self.finished_freq = []
                    self.freq_ant = {}
                    for l in range(2048):self.freq_ant[l] = []
                    break
            #print self.finished_freq 
            if self.finished_freq == []:
                self.finished_freq.append(0)
                self.pkt_len = 0
            #print max(self.finished_freq)    
            return max(self.finished_freq),self.pkt_len     

    def read_addr(self,xeng):
        self.snap_addr[xeng].flush()
        self.snap_addr[xeng].seek(0)
        self.snap_addr[xeng].flush()
        return struct.unpack('L',self.snap_addr[xeng].read(4))[0]

    def get_hw_timestamp(self,xeng):
        #self.vacc_mcnt.flush()
        #self.vacc_mcnt.seek(0)
        #self.vacc_mcnt.flush()
        self.vacc_mcnt_l[xeng].flush()
        self.vacc_mcnt_l[xeng].seek(0)
        self.vacc_mcnt_l[xeng].flush()
        self.vacc_mcnt_h[xeng].flush()
        self.vacc_mcnt_h[xeng].seek(0)
        self.vacc_mcnt_h[xeng].flush()
        #return struct.unpack('>L',self.vacc_mcnt.read(4))[0]
        return struct.unpack('>Q',self.vacc_mcnt_h[xeng].read(4)+self.vacc_mcnt_l[xeng].read(4))[0]

    def read_bram(self,xeng,size):
        """Reads "size" bytes from bram for xengine number xeng"""
        self.snap_bram[xeng].flush()
        self.snap_bram[xeng].seek(0)
        return self.snap_bram[xeng].read(size*4)

    def get_acc_len(self):
        a_l=open('/proc/%i/hw/ioreg/acc_len'%(self.pids[0]),'r')
        acc_len=struct.unpack('L',a_l.read(4))[0]
        a_l.close()
        return acc_len

    def snap_get_new(self,xeng):
        self.snap_en[xeng].seek(0)
        self.snap_en[xeng].flush()
        self.snap_en[xeng].write(struct.pack('L',0))
        self.snap_en[xeng].flush()
        self.snap_en[xeng].seek(0)
        self.snap_en[xeng].flush()
        self.snap_en[xeng].write(struct.pack('L',1))
        self.snap_en[xeng].flush()
        self.snap_en[xeng].seek(0)
        self.snap_en[xeng].flush()
        self.snap_en[xeng].write(struct.pack('L',0))
        self.snap_en[xeng].flush()
        
    def empty_buffers(self):
        print 'Flushing buffers...'
        complete=[]
        total_xeng_vectors=[]

        for xnum in range(self.x_per_fpga):
            complete.append(0)
            total_xeng_vectors.append(0)
            self.snap_get_new(xnum)
            #print 'Requested first snap grab for xeng %i'%xnum

        # Wait for data to become available
        num = 0
        while num == 0:
            time.sleep(.005)
            num = self.read_addr(0)

        while sum(complete) < self.x_per_fpga:
            for xnum in range(self.x_per_fpga):
                time.sleep(.005)
                addr = self.read_addr(xnum)
                #print 'Got addr %i on xeng %i.'%(addr,xnum)
                if addr == 0:
                    complete[xnum]=1
                    #print '\t: %i/%i complete.'%(sum(complete),self.x_per_fpga*len(self.pids)),complete
                else:
                    complete[xnum]=0
                    total_xeng_vectors[xnum] += (addr+1)
                self.snap_get_new(xnum)

        for xnum in range(self.x_per_fpga):
            print '\t: Flushed %i vectors for X engine %i. %i/%i complete.'%(total_xeng_vectors[xnum], self.xeng[xnum], sum(complete),self.x_per_fpga)


    def _tx(self):
        """Continuously transmit correlator data over udp packets."""
        target_pkt_size=(self.payload_len+self.casper_sock.header_size)

        self.empty_buffers()
        self.empty_buffers()
        self.empty_buffers()

        n_integrations = 0

        complete=[]
        timestamp=[]
        rounded_timestamp=[]
        realtime_diff=[]
        realtime_last=[]
        int_xeng_vectors=[]

        for x in range(self.x_per_fpga):
            timestamp.append(self.get_hw_timestamp(x))
            rounded_timestamp.append( (timestamp[x]/self.timestamp_rnd) * self.timestamp_rnd)
            realtime_diff.append(0)
            realtime_last.append(time.time())

        data = []
        for xnum in range(self.x_per_fpga):
            complete.append(0)
            int_xeng_vectors.append(0)
            data.append([])
            self.snap_get_new(xnum)
            #print 'Requested first snap grab for xeng %i'%xnum
        
        #t1 = time.time()
        #self.get_corr_read_missing(self.corr_read_missing)
        #print 'time to get corr_read_missing data and save in redis = ', time.time() - t1
        freq_offset = 0
        self.pkt_len = 0
        while True:
            t_fullspec_start = time.time()
            # Wait for data to become available
            num = 0
            cnt=0
            while num == 0:
                if cnt == 0:
                    freq_offset,self.pkt_len = self.xaui_parse(self.snap_xaui_ram(self.pkt_len,offset=freq_offset))
                    #print 'freq_offset =', freq_offset
                self.adc_amplitudes()
                time.sleep(.1)
                num = self.read_addr(0)
                cnt+= 1
                self.get_corr_read_missing(self.corr_read_missing)
                #print cnt    

            t_fullspec_start_2 = time.time()
            while sum(complete)<self.x_per_fpga:
                for xnum in range(self.x_per_fpga):
                    addr = self.read_addr(xnum)
                    if addr == 0:
                        complete[xnum]=1
                    else:
                        complete[xnum]=0
                        int_xeng_vectors[xnum] += (addr+1)
                        data[xnum].append(self.read_bram(xnum,addr+1))
                    self.snap_get_new(xnum)

            for x in range(self.x_per_fpga):
                timestamp[x] = self.get_hw_timestamp(x)
                realtime_diff[x]=time.time() - realtime_last[x]
                rounded_timestamp[x] = (timestamp[x]/self.timestamp_rnd) * self.timestamp_rnd
                realtime_last[x]=time.time()

            #Now that we have collected all this integration's data, send the packets:
            for xnum in range(self.x_per_fpga):
                print '[%6i] Grabbed %i vectors for X engine %i with timestamp %i (diff %4.2fs).'%(n_integrations,int_xeng_vectors[xnum], self.xeng[xnum], rounded_timestamp[xnum],realtime_diff[xnum])
                data[xnum] = ''.join(data[xnum])

                print 'time to get all data = ', time.time() - t_fullspec_start
                print 'time to get data - sleep =', time.time() - t_fullspec_start_2


                #n_bls=16*17/2
                #bls=2
                #for chan in range(20):
                #    index=(chan*n_bls+bls)*2*4*4
                #    print 'Chan %4i (%4i): '%(chan,index),struct.unpack('>ii',data[xnum][index:index+8])


                for cnt in range((len(data[xnum])/self.payload_len)):
                    if self.casper_sock.send_packet(rounded_timestamp[xnum], self.xeng[xnum], cnt*self.payload_len, 139264, data[xnum][cnt*self.payload_len:(cnt+1)*self.payload_len]) != target_pkt_size: print 'TX fail!' 
                    #time.sleep(0.000001)
                    #print '.',
                print '\n'
                data[xnum]=[]
                complete[xnum]=0
                int_xeng_vectors[xnum]=0
            n_integrations += 1
            if n_integrations == 1:
                try:
                    self.mcache.set('px%d:integration'%(self.xeng[0]+1),0)
                except Exception, e: print 'REDIS ERROR (set integration)', e
                
                
            try:self.mcache.incr('px%d:integration'%(self.xeng[0]+1))
            except Exception, e: print 'REDIS ERROR (integration increment)', e

if __name__ == '__main__':
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('cn_tx.py [options] pid')
    p.set_description(__doc__)
    p.add_option('-i', '--udp_ip', dest='udp_ip', default='192.168.100.1',
        help='IP address to use for UDP transmission of correlator data.  Default is 192.168.100.1')
    p.add_option('-k', '--udp_port', dest='udp_port', type='int', default=7148,
        help='Port to use for UDP correlator data transmission.  Default is 7148')
    p.add_option('-x', '--x_per_fpga', dest='x_per_fpga', type='int', default=2,
        help='Number of X engines per FPGA.  Default is 2')
    p.add_option('-v', '--verbose', dest='verbose', action='store_true',
        help='Be verbose')
    p.add_option('-l', '--payload_len', dest='payload_len', type='int',default=4096,
        help='The length in bytes of each packet (data or payload only). Default 4096')
    p.add_option('-t', '--timestamp_rounding', dest='timestamp_rounding', type='int', default=1024*1024,
        help='Round-off the timestamp to the nearest given value. Default is 1024*1024.')
    p.add_option('-n', '--n_xaui', dest='n_xaui', type=int, default=1,
        help='Number of xaui ports used per fpga.')
    opts, args = p.parse_args(sys.argv[1:])
    if len(args) < 1: 
        print 'Please specify PID of Xengine BORPH process.'
        sys.exit()
    pid =  int(args[0])
    c = CorrTX(pid, ip=opts.udp_ip, x_per_fpga=opts.x_per_fpga, port=opts.udp_port, payload_len=opts.payload_len, timestamp_rnd=opts.timestamp_rounding, verbose=opts.verbose)
    c.start()
