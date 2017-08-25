'''A module implementing high-level interfaces to Casper Correlators'''
import aipy as a, numpy as n
import rx
import os, ephem,time
import ctypes
import hera_mc.utils
import astropy.time
import json

CLOCK_MONOTONIC_RAW = 4 # see <linux/time.h>

class timespec(ctypes.Structure):
    _fields_ = [
        ('tv_sec', ctypes.c_long),
        ('tv_nsec', ctypes.c_long)
    ]

librt = ctypes.CDLL('librt.so', use_errno=True)
clock_gettime = librt.clock_gettime
clock_gettime.argtypes = [ctypes.c_int, ctypes.POINTER(timespec)]

def monotonic_ns():
    t = timespec()
    if clock_gettime(CLOCK_MONOTONIC_RAW, ctypes.pointer(t)) != 0:
        errno_ = ctypes.get_errno()
        raise OSError(errno_, os.strerror(errno_))
    return t.tv_sec * 1e9 + t.tv_nsec

def expand_cminfo(cminfo, field, default='', nants=128):
    x = [default for i in range(nants)]
    for an, ant in enumerate(cminfo['antenna_numbers']):
        x[ant] = cminfo[field][an]
    return x




def start_uv_file(filename, aa, pols, nchan, sfreq, sdf, inttime, cminfo=None):
    # Set the type of 'corr' to 'r' since we asked miriad to store visdata
    # as floats (aka reals) to avoid dynamic range problem with scaled shorts.
    # Note that this requires aipy version 1.0.1 or newer!!!
    uv = a.miriad.UV(filename, status='new', corrmode='r')
    uv._wrhd('obstype','mixed-auto-cross')
    uv._wrhd('history','CORR-DACQ: created file.\n')
    uv.add_var('telescop','a'); uv['telescop'] = 'PAPER'
    uv.add_var('operator','a'); uv['operator'] = 'PAPER'
    uv.add_var('version' ,'a'); uv['version'] = '1.3.0'
    uv.add_var('epoch'   ,'r'); uv['epoch'] = 2000.
    uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
    uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
    uv.add_var('dec'     ,'d'); uv['dec'] = aa.lat
    uv.add_var('obsdec'  ,'d'); uv['obsdec'] = aa.lat
    uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
    uv.add_var('altitude','d'); uv['altitude'] = aa.elev
    uv.add_var('npol'    ,'i'); uv['npol'] = len(pols)
    uv.add_var('nspect'  ,'i'); uv['nspect'] = 1
    uv.add_var('nants'   ,'i'); uv['nants'] = len(aa)
    uv.add_var('antpos'  ,'d')
    antpos = n.array([ant.pos for ant in aa], dtype=n.double)
    uv['antpos'] = antpos.transpose().flatten()
    uv.add_var('sfreq'   ,'d'); uv['sfreq'] = sfreq
    uv.add_var('freq'    ,'d'); uv['freq'] = sfreq
    uv.add_var('restfreq','d'); uv['restfreq'] = sfreq + sdf*nchan/2
    uv.add_var('sdf'     ,'d'); uv['sdf'] = sdf
    uv.add_var('nchan'   ,'i'); uv['nchan'] = nchan
    uv.add_var('nschan'  ,'i'); uv['nschan'] = nchan
    uv.add_var('inttime' ,'r'); uv['inttime'] = float(inttime)
    if cminfo is not None:
        uv.add_var('antnames', 'a'); uv['antnames'] = '[' + ', '.join(expand_cminfo(cminfo, 'antenna_names')) + ']'
        uv.add_var('cmver',    'a'); uv['cmver']    = str(cminfo['cm_version'])
        uv.add_var('st_type',  'a'); uv['st_type']  = '[' + ', '.join(expand_cminfo(cminfo, 'station_types')) + ']'
        try:
            cminfo.pop('antenna_positions') # this np array can't be json serialized
            uv.add_var('cminfo', 'a'); uv['cminfo']   = json.dumps(cminfo)
        except:
            print 'failed to jsonify cminfo'
    # These variables just set to dummy values
    uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
    uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
    uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
    uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.
    # These variables will get updated every spectrum
    uv.add_var('coord'   ,'d')
    uv.add_var('time'    ,'d')
    uv.add_var('lst'     ,'d')
    uv.add_var('ra'      ,'d')
    uv.add_var('obsra'   ,'d')
    uv.add_var('baseline','r')
    uv.add_var('pol'     ,'i')
    # Variables added in July 2017 for HERA deployment
    uv.add_var('startt'   ,'d') # double precision gps time of first integration
    uv.add_var('stopt'    ,'d') # double precicion gps time of last integration
    uv.add_var('obsid'    ,'i') # int(floor(starttime))
    uv.add_var('duration' ,'i') # integer seconds of an observation (starttime - stoptime)
    return uv

class DataReceiver(rx.BufferSocket):
    def __init__(self, aa, nants_per_feng=4, pols=['xx','yy','xy','yx'], adc_rate=100000000.,
            nchan=2048, xeng_chan_mode=0, sfreq=0.121142578125, sdf=7.32421875e-05,
            inttime=14.3165578842, t_per_file=ephem.hour, payload_data_type=0,
            nwin=4, bufferslots=128, payload_len=8192, sdisp=0, sdisp_destination_ip="127.0.0.1",
            acc_len=1024*128, redis=None, hookup=None, cminfo=None):
        rx.BufferSocket.__init__(self, item_count=bufferslots, payload_len=payload_len)
        self.cb = rx.CollateBuffer(nant=len(aa), nants_per_feng=nants_per_feng, npol=len(pols),
            nchan=nchan, xeng_chan_mode=xeng_chan_mode, nwin=nwin, sdisp=sdisp,
            sdisp_destination_ip=sdisp_destination_ip, acc_len=acc_len)
        # Define a file-writing callback that starts/ends files when
        # appropriate and updates variables
        npol = len(pols)
        self.uv = None
        self.adc_rate=float(adc_rate)
        self.filestart = 0.
        self.current_time = 0
        self.t_per_file = t_per_file
        self.adc_rate=float(adc_rate)
        self.redis = redis
        self.hookup = hookup
        self.uvio_total_ns = 0
        self.uvio_calls = 0
        self.cminfo = cminfo

        def filewrite_callback(i,j,pol,tcnt,data,flags):
            # The parameters i, j, and pol have been derived solely from
            # correlator input index and may have nothing to do with actual
            # antpol information.  If we have hookup information, we can derive
            # the correlator input indexs for both inputs and use that to
            # lookup the two antpols in this baseline and relabel it
            # accordingly.  The correlator has conjugated the data according to
            # (i,pol_0) < (j,pol_1), but the data may need to be conjgated
            # (again) depending on how i, j, and pol map to reality.  See the
            # comments below for more details.
            if self.hookup is not None:
                # If we get any exceptions, skip the hookup stuff.  This is
                # likely to lead to crap data, but it's better than crashing
                # and leading to no data.
                try:
                    # Compute input indexes
                    in0 = 2*i + (ord(pols[pol][0]) - ord('x'))%2
                    in1 = 2*j + (ord(pols[pol][1]) - ord('x'))%2
                    # Lookup antnums and pol strings.  The hookup dictionary
                    # maps intput number to [antnum, p], where antnum is
                    # integer and p is lowercase 'x' or 'y'.
                    a0, p0 = hookup[in0]
                    a1, p1 = hookup[in1]
                    # Lookup pol-pair index ("p" is "reality" polarization
                    # index, "pol" is "correlator computed" polarization.
                    p = pols.index(p0 + p1)

                    # Conjugate the data if needed.  The following chart shows
                    # when conjugation is needed (CONJ) and when it's not (OK):
                    #
                    #                   +---------------------+
                    #                   | Correlator Computed |
                    #                   +------+------+-------+
                    #                   | Auto | Auto | Cross |
                    #                   |  XY  |  YX  |  AB   |
                    #     +---+---------+------+------+-------+
                    #     | R | Auto XY |  OK  | CONJ |  OK   | Row1
                    #     | e +---------+------+------+-------+
                    #     | a | Auto YX | CONJ |  OK  | CONJ  | Row2
                    #     | l +---------+------+------+-------+
                    #     | i |Cross CD |  OK  | CONJ |  OK   | Row3
                    #     | t +---------+------+------+-------+
                    #     | y |Cross DC | CONJ |  OK  | CONJ  | Row4
                    #     +---+---------+------+------+-------+
                    #                     Col1   Col2   Col3
                    #
                    #  Notice that Col2 is the inverse of Col1 and Col3.

                    need_conj = False

                    # If Reality "Cross CD" or "Auto YX" (i.e. Row4 or Row2)
                    if (a0 > a1) or (a0 == a1 and pols[p] == 'yx'):
                        need_conj = not(need_conj)

                    # If Correlator Computed "Auto YX" (i.e. Col2)
                    if i == j and pols[pol] == 'yx':
                        need_conj = not(need_conj)

                    if need_conj:
                        # Do in-place conjugation
                        data.conj(data)

                    # Ensure that the baseline is ordered properly in MIRIAD
                    if a0 > a1:
                        # Swap antennas and pols
                        a0, a1 = a1, a0
                        p = pols.index(p1 + p0)

                    # Update i, j, and pol with new values
                    i, j, pol = a0, a1, p
                except:
                   print sys.exc_info()[0]
                   print 'Ignoring hookup info!!!'
                   self.hookup = None

            # Update time and baseline calculations if tcnt changes, possibly
            # ending a file and starting a new one if necessary

            #t is julian date. ephem measures days since noon, 31st Dec 1899 (!random!), so we need an offset from uni time:
            t_unix = tcnt / self.adc_rate
            t = a.phs.ephem2juldate((t_unix + 2209032000)*ephem.second)
            

            #print "fwrite callback:",i,j,pol,tcnt,data.size,flags.size

            #if i==0 and j==0 and pol==0: print '0-0-0: ',sum(data)
            #if i==0 and j==0 and pol==1: print '0-0-1: ',sum(data)
            #if i==1 and j==1 and pol==0: print '1-1-0: ',sum(data)
            #if i==1 and j==1 and pol==1: print '1-1-1: ',sum(data)

            if (t != self.current_time):
                self.current_time = t

                if self.uvio_calls > 0 and pol == 0:
                    print 'Total time spent in %s uvio calls is %.1f seconds (%f ns/call)' % (
                        self.uvio_calls, float(self.uvio_total_ns)/1e9, float(self.uvio_total_ns)/self.uvio_calls)
                    self.uvio_calls = 0
                    self.uvio_total_ns = 0

                if (t > (self.filestart + self.t_per_file)) or self.uv == None:
                    # Get the astropy.Time version of the time -- we're going to need it
                    # to write the time-related file variables
                    astro_time = astropy.time.Time(t_unix, format='unix')
                    if self.uv != None:
                        # Write the stop-time variables
                        uv['stoptime'] = astro_time.gps
                        uv['duration'] = int(n.floor(astro_time.gps - uv['starttime']))
                        # Work with filename for the given polarization
                        # Get reference to the UV object self.uv and set
                        # self.uv to None so that the UV object can be
                        # deleted without altering the ordering of other
                        # elements in self.uv array.  This is insane but
                        # seemingly necessary since the UV destructor *MUST* be
                        # called to properly close the dataset, but if one were
                        # to do the "normal" thing of "del(self.uv)", then
                        # the elements in self.uv after pol are shifted down
                        # and things get all confused.
                        uv, self.uv = self.uv, None
                        del(uv)
                        print 'Ending file:',
                        print self.fname, '->', self.fname.replace('.tmp','')
                        os.rename(self.fname, self.fname.replace('.tmp',''))
                    self.fname = 'zen.%07.5f.uv.tmp' % t
                    self.filestart = t
                    print a.phs.juldate2ephem(t),
                    print 'Starting file:', self.fname
                    self.uv = start_uv_file(
                        self.fname, aa, pols=pols, nchan=nchan,
                        sfreq=sfreq, sdf=sdf, inttime=inttime, cminfo=self.cminfo)
                    # Set the file start time and obsid
                    self.uv['obsid'] = hera_mc.utils.calculate_obsid(astro_time)
                    self.uv['startt'] = astro_time.gps

                aa.set_jultime(t)
                lst = aa.sidereal_time()
                self.uv['lst'] = lst
                self.uv['ra'] = lst
                self.uv['obsra'] = lst

            self.uv['pol'] = a.miriad.str2pol[pols[pol]]
            crd = aa[j].pos - aa[i].pos
            preamble = (crd, t, (i,j))

            # Only clip RFI if visibilities are being stored as scaled shorts
            # and it is not an autocorrelation.
            if self.uv.vartable['corr'] == 'j' and (i!=j or pols[pol] in ['xy','yx']):
                dabs = n.abs(data)
                # Clips RFI to amplitude 1, motivated by desire to avoid miriad
                # dynamic range readout issue for scaled shorts.
                data = n.where(dabs>1.0,data/dabs,data)

            # If redis is given and this is an xx or yy autocorrelation
            if self.redis is not None and i==j and pols[pol] in ['xx','yy']:
                try:
                    key = 'visdata://%d/%d/%s' % (i, j, pols[pol])
                    self.redis.hmset(key,
                        {'time': str(t), 'data': n.abs(data).tostring()})
                except Exception, e:
                    print 'redis exception: %s' % e

            self.uvio_total_ns -= monotonic_ns()
            self.uv.write(preamble, data, flags)
            self.uvio_total_ns += monotonic_ns()
            self.uvio_calls += 1

        self.cb.set_callback(filewrite_callback)
        self.set_callback(self.cb)

    def set_start_time(self, start_jd, tscale):
        self.start_jd = start_jd
        self.tscale = tscale

    def stop(self):
        rx.BufferSocket.stop(self)
        # Close data sets
        if self.uv != None:
            # Delete entire self.uv array, which will delete each element,
            # which will call their destructors (hopefully, destructors are not
            # guaranteed to be called!), which will close the datasets.  In
            # practice this works OK with the standard python implementation
            # (aka "CPython").
            del(self.uv)
        # Rename datasets
        if self.fname is not None:
            print 'Ending file:',
            print self.fname, '->', self.fname.replace('.tmp','')
            os.rename(self.fname, self.fname.replace('.tmp',''))
