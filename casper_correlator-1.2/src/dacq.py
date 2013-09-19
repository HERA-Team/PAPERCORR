'''A module implementing high-level interfaces to Casper Correlators'''
import aipy as a, numpy as n
import rx
import os, ephem,time

def start_uv_file(filename, aa, pols, nchan, sfreq, sdf, inttime):
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
    return uv

class DataReceiver(rx.BufferSocket):
    def __init__(self, aa, nants_per_feng=4, pols=['xx','yy','xy','yx'], adc_rate=100000000.,
            nchan=2048, xeng_chan_mode=0, sfreq=0.121142578125, sdf=7.32421875e-05,
            inttime=14.3165578842, t_per_file=ephem.hour, payload_data_type=0,
            nwin=4, bufferslots=128, payload_len=8192, sdisp=0, sdisp_destination_ip="127.0.0.1",
            acc_len=1024*128, redis=None):
        rx.BufferSocket.__init__(self, item_count=bufferslots, payload_len=payload_len)
        self.cb = rx.CollateBuffer(nant=len(aa), nants_per_feng=nants_per_feng, npol=len(pols),
            nchan=nchan, xeng_chan_mode=xeng_chan_mode, nwin=nwin, sdisp=sdisp,
            sdisp_destination_ip=sdisp_destination_ip, acc_len=acc_len)
        # Define a file-writing callback that starts/ends files when
        # appropriate and updates variables
        self.uv = None
        self.adc_rate=float(adc_rate)
        self.filestart = 0.
        self.current_time = 0
        self.t_per_file = t_per_file
        self.redis = redis

        def filewrite_callback(i,j,pol,tcnt,data,flags):
            # Update time and baseline calculations if tcnt changes, possibly
            # ending a file and starting a new one if necessary

            #t is julian date. ephem measures days since noon, 31st Dec 1899 (!random!), so we need an offset from uni time:
            t = a.phs.ephem2juldate(((tcnt/self.adc_rate) + 2209032000)*ephem.second)

            #print "fwrite callback:",i,j,pol,tcnt,data.size,flags.size

            #if i==0 and j==0 and pol==0: print '0-0-0: ',sum(data)
            #if i==0 and j==0 and pol==1: print '0-0-1: ',sum(data)
            #if i==1 and j==1 and pol==0: print '1-1-0: ',sum(data)
            #if i==1 and j==1 and pol==1: print '1-1-1: ',sum(data)

            if (t != self.current_time):
                self.current_time = t

                if (t > (self.filestart + self.t_per_file)) or self.uv == None:
                    if self.uv != None:
                        del(self.uv)
                        print 'Ending file:',
                        print self.fname, '->', self.fname.replace('.tmp','')
                        os.rename(self.fname, self.fname.replace('.tmp',''))
                    self.fname = 'zen.%07.5f.uv.tmp' % t
                    print a.phs.juldate2ephem(t),
                    print 'Starting file:', self.fname
                    self.uv = start_uv_file(self.fname, aa, pols=pols,
                        nchan=nchan, sfreq=sfreq, sdf=sdf, inttime=inttime)
                    self.filestart = t

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
                    self.redis.hmset(key, {'time': str(t), 'data': n.abs(data).tostring()})
                except Exception, e:
                    print 'redis exception: %s' % e

            self.uv.write(preamble, data, flags)

        self.cb.set_callback(filewrite_callback)
        self.set_callback(self.cb)

    def set_start_time(self, start_jd, tscale):
        self.start_jd = start_jd
        self.tscale = tscale
    def stop(self):
        rx.BufferSocket.stop(self)
        if self.uv != None:
            del(self.uv)
            print 'Ending file:',
            print self.fname, '->', self.fname.replace('.tmp','')
            os.rename(self.fname, self.fname.replace('.tmp',''))



