
scatPath = 'D:\LISST-VSF\Measurements\Line P Feb 2018\V0540037.VSF';
%zscPath = 'D:\LISST-VSF\Measurements\Line P June 2017\Z1571050.VSF';
zscPath = 'D:\LISST-VSF\FactoryBackground_SN1662.VSF';
alpha = [];
plotData = 1; %0 off or 1 on for plotting
SN = 1662;
proc = lisstvsf_makep(SN, scatPath, zscPath, alpha, plotData);