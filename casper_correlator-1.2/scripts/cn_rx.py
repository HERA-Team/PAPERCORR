#!/usr/bin/python -u
import casper_correlator,corr,ephem,aipy,numpy,sys,socket,time,struct,syslog,signal

syslog.openlog('cn_rx.py')

# 2-14-2011 Z.A. added 16-31 in 'ants'. preparation for 64 input corr.
if sys.argv[1:]==[]:
    print 'Please specify configuration file.'
    exit()      
lh=corr.log_handlers.DebugLogHandler()
c=corr.corr_functions.Correlator(sys.argv[1],lh)
nants = c.config['n_ants']
nants_per_feng = c.config['n_ants_per_feng']
port = c.config['rx_udp_port']
n_chans = c.config['n_chans']
xeng_chan_mode = c.config['xeng_chan_mode']
bandwidth = c.config['adc_clk']/2 # GHz
sdf = bandwidth/n_chans
sfreq = bandwidth # Second Nyquist zone
location=0,0,0
acc_len = c.config['acc_len'] * c.config['xeng_acc_len']
int_time = 2*n_chans*acc_len/(bandwidth*2*1e9) #integration time in seconds
#acc_len = 1 
 # incoming data divided by this number for correct scaling
t_per_file=ephem.minute*10
n_windows_to_buffer=4
n_bufferslots=10240
max_payload_len=8192
payload_data_type = c.config['payload_data_type']
ants=[]
for i in range(nants):
    ants.append((i,i,i))
pols=['xx','yy','xy','yx']
freqs = numpy.arange(n_chans, dtype=numpy.float) * sdf + sfreq
beam = aipy.phs.Beam(freqs)
ants = [aipy.phs.Antenna(a[0],a[1],a[2],beam) for a in ants]
aa = aipy.phs.AntennaArray(ants=ants, location=location)
sdisp_destination_ip = "127.0.0.1"
rx=None
use_redis = c.redis

if len(sys.argv) > 2 and sys.argv[2] == '--no-redis':
    use_redis = None

if use_redis is None:
    print 'Autocorrelations will NOT be stored in redis'
else:
    print 'Autocorrelations will be stored in redis'

# Function to handle SIGINT
def stop_taking_data(*args):
    global rx
    if rx is not None:
      print 'stopping cn_rx.py'
      rx.stop()
    exit(0)

try:
    rx=casper_correlator.dacq.DataReceiver(aa, nants_per_feng=nants_per_feng,
                pols=pols, adc_rate=100000000, nchan=n_chans,
                xeng_chan_mode=xeng_chan_mode, sfreq=sfreq, sdf=sdf,
                inttime=int_time, t_per_file=t_per_file,
                nwin=n_windows_to_buffer, bufferslots=n_bufferslots, 
                payload_len=max_payload_len, payload_data_type=payload_data_type,
                sdisp=0, sdisp_destination_ip=sdisp_destination_ip,
                acc_len=acc_len, redis=use_redis)
    rx.start(port)

    signal.signal(signal.SIGINT, stop_taking_data)
    
    time.sleep(5)

    print 'Setting time lock...'    

    # Try new name first, then old name
    try:
        trig_time = float(c.redis.get('roachf_init_time'))
    except:
        trig_time = float(c.redis.get('mcount_initialize_time'))
    time_skt = socket.socket(type=socket.SOCK_DGRAM)
    pkt_str=struct.pack('>HHHHQ',0x5453,3,0,1,int(trig_time))
    time_skt.sendto(pkt_str,(c.config['rx_udp_ip_str'],c.config['rx_udp_port']))
    time_skt.close()
    print 'Time Pkt sent...'

    while True:
      #capture a bunch of stuff here
      time.sleep(10)

except(KeyboardInterrupt):
    stop_taking_data()
