main code: get_lisst_vsf.v4.m
The input of beamc is optional. If you have independent measurements of the beam attenuation coefficient, I suggest to use it as an input. If no, just set beamc as empty. 
The output of the main code is a structure variable 'proc'.
proc.Angles: scattering angles of the LISST-VSF. The first 32 angles are for the LISST unit and the left 141 angles are for the eyeball component.
proc.pmt_gain: PMT voltage
proc.P11: VSFs of both the LISST unit and the eyeball component. The VSFs of the eyeball component are processed with default method.
proc.P11_absCal: VSFs of both the LISST unit and the eyeball component. The VSFs of the eyeball component are processed with the absolute method (Hu et al, 2019 OE).
proc.P11_opt: VSFs of both the LISST unit and the eyeball component. The VSFs of the eyeball component are processed with the default method if PMT voltage<=528mV or they are processed with the absolute method. 