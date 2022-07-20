

function calfact = lisstvsf_getcal(SN)

switch SN  
    case 1425
        
        calfact.SN = 1425;
        calfact.depth_slope = 0.01;
        calfact.depth_offset = 0;
        calfact.temp_slope = 0.01;
        calfact.temp_offset = 0;
        calfact.bat_slope = 0.01;
        calfact.bat_offset = 0;
        calfact.angle_offset = -2;
        calfact.Watt_per_count_on_rings = 1.9e-10; 
        calfact.Watt_per_count_laser_ref = 3.76e-6 .* 1.04; % multiplied by 1.04 for Fresnel loss
        calfact.ND = 10.23; 
        calfact.laser_power_change_factor = 9; 
        calfact.HWPlate_transmission = 0.99;
        calfact.focal_length = 53;
       % calfact.geometric_cal_coeff = [5.74393643093325e-09,-2.47087923775527e-06,0.000377167939477837,-0.0255620583681450,1.44582831425705];
       calfact.geometric_cal_coeff = [0 1];
        
        % Geometric dcal
        calfact.dcal = [1.01790833330000	0.992134894300000	1.01081606840000	0.994928828300000	1.00437070280000	0.998918404700000	0.998590554600000	1.00420485260000	1.00007630710000	0.998899970300000	1.00094967890000	1.00040191990000	1.00111296400000	1.00046768960000	1.02135541230000	0.999902616600000	1.11156301020000	1.12066675740000	1.24936993590000	1.16431994250000	1.33556567670000	1.20908922660000	1.07815400160000	1.76207523180000	1.55085634300000	2.53041187530000	2.56385922160000	3.57572117650000	3.96319866420000	5.01664107560000	5.77011140700000	7.01385046680000];
        
        % De-vignetting factors
        calfact.dvig = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1.02101205238326,1.20581086480367,1.43697017345649,1.69141875258289,2.00024365577338,2.44936150814052];
         
    otherwise
        calfact = [];
        warning(['No calibration factors for instrument ' num2str(SN) '.']);
        
end

