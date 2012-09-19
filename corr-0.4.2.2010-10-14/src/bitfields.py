import construct, corr_functions

register_fengine_control = construct.BitStruct('feng_ctl',
        construct.Flag('gbe_gpu_rst'),        #31
        construct.Flag('gbe_sw_rst'),         #30
        construct.Flag('loopback_mux_rst'),   #29
        construct.Flag('cnt_rst'),            #28
        construct.BitField('fft_preshift',2), #26-27
        construct.BitField('gpio_monsel',3),  #23,24,25
        construct.Flag('fft_tvg2'),           #22
        construct.Flag('fft_tvg1'),           #21
        construct.Flag('gbe_gpu_disable'),    #20
        construct.Flag('use_qdr_tvg'),        #19
        construct.Flag('gbe_sw_disable'),     #18
        construct.Flag('arm_rst'),            #17
        construct.Flag('sync_rst'),           #16
        construct.Flag('arm_noise'),          #15
        construct.Padding(12),                #3-14
        construct.Flag('lb_err_cnt_rst'),     #2
        construct.Padding(2))                 #0-1

