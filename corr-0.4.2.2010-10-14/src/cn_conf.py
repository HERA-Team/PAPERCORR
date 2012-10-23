import ConfigParser,exceptions,socket,struct
"""
Library for parsing CASPER correlator configuration files

Author: Jason Manley
Revision: 2008-02-08

Revs:
2009-12-10: JRM Changed IP address formats to input strings, but return integers.
                Added __setitem__, though it's volatile.
                added calculation of n_bls, bandwidth, integration time etc.
2008-02-08: JRM Replaced custom tokeniser with string.split
                Changed EQ to EQ_polys
                Changed max_payload_len to rx_buffer_size

"""
LISTDELIMIT = ','
PORTDELIMIT = ':'


class CorrConf:    
    def __init__(self, config_file):
        self.config_file=config_file
#        self.cfp=open(config_file,'rwb')
        self.cfp=open(config_file,'rb')
        self.cp = ConfigParser.ConfigParser()
        self.cp.readfp(self.cfp)
        self.config=dict()
        self.read_all()

    def __getitem__(self,item):
        return self.config[item]

    def __setitem__(self,item,value):
        self.config[item]=value

    def file_exists(self):
        try:
            f = open(self.config_file)
        except IOError:
            exists = False
        else:
            exists = True
            f.close()
        return exists

    def read_all(self):
        if not self.file_exists():
            raise RuntimeError('Error opening config file')
        self.config['correlator']=dict()
    
        #Get the servers:
        self.config['servers']=[]
        serverstr= self.cp.get('borphserver','servers')

        #print 'got serverstring ',serverstr
        server_port_strings = serverstr.split(LISTDELIMIT)
        self.servers=[]
        for i in range(len(server_port_strings)):
            tokens_out = server_port_strings[i].split(PORTDELIMIT)
            if len(tokens_out)==2:
                self.config['servers'].append({'server':tokens_out[0]})                
                self.config['servers'][i].update({'port':int(tokens_out[1])})
            else:
                raise RuntimeError('Error parsing server section: did not get port and ip addr for server') 

        self.config['timeserver'] = self.cp.get('borphserver','timeserver')
        self.config['bitstream'] = self.cp.get('borphserver','bitstream')

        #get the correlator config stuff:
        self.read_int('correlator','xmode')
        self.read_int('correlator','n_chans')
        self.read_int('correlator','n_ants')
        self.read_int('correlator','n_ants_per_feng')
        self.read_int('correlator','fft_shift')
        self.read_str('correlator','trig_mode')
        self.read_str('correlator','int_trig_serial_port')
        self.read_int('correlator','acc_len')
        self.read_float('correlator','adc_clk')
        self.read_int('correlator','n_stokes')
        self.read_int('correlator','x_per_fpga')
        self.read_int('correlator','n_ants_per_xaui')
        self.read_int('correlator','clk_per_sync')
        self.read_int('correlator','xeng_acc_len')
        self.read_int('correlator','xeng_chan_mode')
        self.read_float('correlator','ddc_mix_freq')
        self.read_int('correlator','ddc_decimation')
        self.read_int('correlator','10gbe_port')
        self.read_int('correlator','10gbe_sw_port')
        self.read_int('correlator','10gbe_gpu_port')
        self.read_int('correlator','10gbe_pkt_len')
        self.read_int('correlator','feng_bits')
        self.read_int('correlator','feng_fix_pnt_pos')
        self.read_int('correlator','x_eng_clk')
        self.read_int('correlator','n_xaui_ports_per_fpga')
        self.read_int('correlator','adc_bits')
        self.read_int('correlator','adc_levels_acc_len')
        self.read_int('correlator','antenna_offset_addr')
        self.read_int('correlator','payload_data_type')
        self.config['10gbe_ip']=struct.unpack('>I',socket.inet_aton(self.cp.get('correlator','10gbe_ip')))[0]
        self.config['10gbe_sw_ip']=struct.unpack('>I',socket.inet_aton(self.cp.get('correlator','10gbe_sw_ip')))[0]
        self.config['10gbe_gpu_ip']=struct.unpack('>I',socket.inet_aton(self.cp.get('correlator','10gbe_gpu_ip')))[0]
        #print '10GbE IP address is %i'%self.config['10gbe_ip']

        gpuips = self.cp.get('correlator','gpu_ip').split(LISTDELIMIT)
        if len(gpuips) != len(self.config['servers']):
            print 'Number of fengines (servers) and xengines (gpu_ip) does not match!'
            self.config['gpu_ips'] = []
        else:
            self.config['gpu_ips'] = []
            for ip in gpuips:
                self.config['gpu_ips'].append(struct.unpack('>I', socket.inet_aton(socket.gethostbyname(ip)))[0])

        self.config['n_bls']=self.config['n_ants']*(self.config['n_ants']+1)/2
        if self.config['ddc_mix_freq']<=0:
            #We're dealing with a "real" PFB, either wideband or narrowband.
            self.config['bandwidth']=self.config['adc_clk']*1000000000./2
            self.config['center_freq']=self.config['adc_clk']*1000000000./4
            self.config['int_time']= float(self.config['n_chans'])*self.config['xeng_acc_len']*self.config['acc_len']/(self.config['bandwidth'])
        else:
            #We're dealing with a complex PFB with a DDC upfront.
            self.config['bandwidth']=self.config['adc_clk']*1000000000./self.config['ddc_decimation']
            self.config['center_freq']=self.config['adc_clk']*1000000000.*self.config['ddc_mix_freq']
            self.config['int_time']= float(self.config['n_chans'])*self.config['xeng_acc_len']*self.config['acc_len']/(self.config['adc_clk']*1000000000./self.config['ddc_decimation'])

        self.config['sync_time']=-1

        #get the receiver section:
        self.config['receiver']=dict()
        self.read_int('receiver','rx_udp_port')
        self.read_int('receiver','rx_pkt_payload_len')
        self.read_int('receiver','instance_id')
        self.config['rx_udp_ip_str']=self.cp.get('receiver','rx_udp_ip')
        self.config['rx_udp_ip']=struct.unpack('>I',socket.inet_aton(self.cp.get('receiver','rx_udp_ip')))[0]
        #print 'RX UDP IP address is %i'%self.config['rx_udp_ip']

        #equalisation section:
        self.read_bool('equalisation','auto_eq')
        self.read_int('equalisation','eq_decimation')
        #self.read_int('equalisation','eq_brams_per_pol_interleave')
        self.config['eq']=dict()
        self.config['eq']['eq_polys']=[]
        if not self.config['xmode']:
            for ant in range(self.config['n_ants']):
                for pol in ['x','y']:
                    ant_eq_str=self.get_line('equalisation','eq_poly_%i%c'%(ant,pol))
                    if (ant_eq_str):
                        self.config['eq']['eq_poly_%i%c'%(ant,pol)]=[int(coef) for coef in ant_eq_str.split(LISTDELIMIT)]
                    else:
                        raise('ERR eq_poly_%i%c'%(ant,pol))

#        # get the antenna info:
#        self.config['antennas']=dict()
#        self.config['antennas']['pos']=[]
#        for ant in range(self.config['n_ants']):
#            ant_pos_str=self.get_line('antennas','pos_%i'%(ant))
#            if (ant_pos_str):
#                self.config['antennas']['pos'].append([float(coef) for coef in ant_pos_str.split(LISTDELIMIT)])
#            else:
#                raise RuntimeError('ERR pos_%i'%(ant))
#
#        self.config['antennas']['location']=[]
#        ant_loc_str=self.get_line('antennas','location')
#        if (ant_loc_str):
#            location = ant_loc_str.split(LISTDELIMIT)
#            if len(location) == 3:
#                location[2]=float(location[2])
#                self.config['antennas']['location']= location
#            else: 
#                raise RuntimeError('ERR setting array location')
#        else:
#            raise RuntimeError('ERR setting array location')

#        self.config['antennas']['order']=[]
#        ant_ord_str=self.get_line('antennas','order')
#        if (ant_ord_str):
#            order = ant_ord_str.split(LISTDELIMIT)
#            if len(order) == self.config['n_ants']:
#                self.config['antennas']['order']= [int(ant) for ant in order]
#            else: 
#                raise RuntimeError('ERR Retrieving antenna ordering')
#        else:
#            raise RuntimeError('ERR Retrieving antenna ordering')

#    def write(self,section,variable,value):
#        self.cp.set(section,variable,str(value))
#        self.cfp.seek(0)
#        self.cp.write(self.cfp)
#        cfp.flush()

    def get_line(self,section,variable):
        return self.cp.get(section,variable)

    def read_int(self,section,variable):
        self.config[variable]=self.cp.getint(section,variable)

    def read_bool(self,section,variable):
        self.config[variable]=(self.cp.get(section,variable) != '0')

    def read_str(self,section,variable):
        self.config[variable]=self.cp.get(section,variable)

    def read_float(self,section,variable):
        self.config[variable]=self.cp.getfloat(section,variable)
