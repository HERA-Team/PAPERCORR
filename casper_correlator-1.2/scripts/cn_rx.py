#!/usr/bin/python -u
import casper_correlator,corr,ephem,aipy,numpy,sys,socket,time,struct,syslog,signal
import yaml

syslog.openlog('cn_rx.py')

# 2-14-2011 Z.A. added 16-31 in 'ants'. preparation for 64 input corr.
if len(sys.argv) < 2:
    print 'Please specify configuration file.'
    exit()

# Returns hookup dictionary mapping correlator input number to [antnum, pol],
# where antnum ranges from 0 to nants-1 and pol is either 'X' or 'Y'.
# Currently gets the data from /etc/papercfg/psa128/fxin_to_antpol.yml.
def get_hookup():
    hookup_name = '/etc/papercfg/psa128/fxin_to_antpol.yml'
    hookup = {}

    try:
        print 'Using fxin-to-antpol file', hookup_name

        hookup_file = open(hookup_name)
        try:
            yaml_hookup = yaml.load(hookup_file.read())
        finally:
            hookup_file.close()

        for k,v in yaml_hookup.items():
            if v is not None:
                # k is "fFAC", where F is '1'..'8', A is 'A'..'H', and C is
                # '1'..'4'.
                #
                # v is "aNNNP", where NNN is a one to three digit number
                # and P is 'X' or 'Y'.

                # Convert k to fx input number.
                f = int(k[1]) - 1
                a = ord(k[2]) - ord('A')
                c = int(k[3]) - 1
                fxin = 32*f + 4*a + c

                # Convert v to [integer_ant, lowecase_pol]
                antpol = [int(v[1:-1]), str.lower(v[-1])]

                # Store in hookup dictionary
                hookup[fxin] = antpol
    except:
        print sys.exc_info()[0]
        print 'Defaulting hookup info!!!'
        hookup = None

    return hookup

# Returns list of 3-tuples of antenna positions for ants 0 to nants-1.
# Currently gets the data from /etc/papercfg/psa128/ant_to_pos.yml.
#
# This function assumes the antenna positions in the YAML file are in [E,N,U]
# (topocentric), but it returns the antenna positions in [N,E,U] (topocentric)
# format because that is just a rotation around the Y/E axis (by latitude) away
# from being in [X,Y,Z] like MIRIAD expects.  Rotation around the Y/E axis by
# latitude is NOT yet done by this function.
#
# This function returns the same units as are in the YAML file.

# This function subtracts the array reference position from the antenna
# positions, so the antenna positions in MIRIAD will all be relative to the
# array reference position.
#
def get_antpos(nants):

    antpos_name = '/etc/papercfg/psa128/ant_to_pos.yml'
    antpos=[]
    # Default to (0,0,0)
    for i in range(nants):
        antpos.append((0,0,0))

    try:
        print 'Using ant-to-pos file', antpos_name
        antpos_file = open(antpos_name)
        try:
            yaml_antpos = yaml.load(antpos_file.read())
        finally:
            antpos_file.close()

        # Setup ref positions (antenna -1) to subtract from positions
        ref = [0, 0, 0]
        if 'a-1' in yaml_antpos:
            ref = yaml_antpos['a-1']

        for i in range(nants):
            # Lookup position for antenna i
            key = 'a%d' % i
            if key in yaml_antpos:
                # Convert from ENU to NEU and store in antpos list as tuple.
                pos = yaml_antpos['a%d'%i]
                antpos[i] = (pos[1]-ref[1], pos[0]-ref[0], pos[2]-ref[2])
            else:
                antpos[i] = (0, 0, 0)
    except:
        print sys.exc_info()[0]
        print 'Defaulting antpos info!!!'
        # Default to (i,i,i) (i.e. the old behavior)
        for i in range(nants):
            antpos.append((i,i,i))

    return antpos

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
ants=get_antpos(nants)
pols=['xx','yy','xy','yx']
freqs = numpy.arange(n_chans, dtype=numpy.float) * sdf + sfreq
beam = aipy.phs.Beam(freqs)
ants = [aipy.phs.Antenna(a[0],a[1],a[2],beam) for a in ants]
aa = aipy.phs.AntennaArray(ants=ants, location=location)
sdisp_destination_ip = "127.0.0.1"
rx=None
use_redis = c.redis
hookup = get_hookup()

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
                acc_len=acc_len, redis=use_redis, hookup=hookup)
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
