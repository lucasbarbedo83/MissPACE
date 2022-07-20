function proc = get_lisst_vsf(scatfile,fbfile,isplot)
%%%

SN = 1662; alpha = []; plotData = 0; %0 off or 1 on for plotting
% absCal = 9.3532e-5; % calibrate p11_uncal to VSF (m-1)
% absCal = 8.823e-5; %  new calibrate coef at 15 degree
absCal_file='D:\SWScatering\LISST-VSF\Measurements\results\Calibration\absCal.txt';
absCal = (dlmread(absCal_file))';

pmt_gain_set = [400 435 470 510 550 595 645 700 760 820];
proc = mylisstvsf_makep(SN, scatfile, fbfile, alpha, plotData);          
gain_value = median(proc.pmt_gain);
 
 idx = gain_value==proc.pmt_gain;
 amp = find(gain_value==pmt_gain_set);
 p11_uncal = proc.p11_uncal(idx,1:136);
 p11_raw = proc.p11_raw(idx,1:136);
 
 nrow = sum(idx);
 absCal = repmat(absCal,nrow,1);
 p11_absCal = p11_uncal./2^(amp-1).*absCal; % covert p11 counts to VSF without water contribution
 p11_absCal_withbgd = p11_raw./2^(amp-1).*absCal; % covert p11 counts to VSF including water contribution
%  p11_absCal_withbgd = p11_raw.*(400/gain_value)^8.6.*absCal;
 
 P11_lisst = nanmean(proc.P11);
 p12_lisst = nanmean(proc.p12); 
  
 proc.p11_absCal = p11_absCal;
 proc.p11_absCal_withbgd = p11_absCal_withbgd;
 proc.P11_absCal = [P11_lisst(1:32) mean(p11_absCal)];
 proc.P11_absCal_withbgd = [P11_lisst(1:32) mean(p11_absCal_withbgd)];
 
 proc.P11_lisst = P11_lisst;
 proc.p12_lisst = p12_lisst;
 
 if isplot
    fh=figure('visible','on'); 
    set(fh, 'Position', [0 0 700 700]);  
    angle = proc.Angles(33:end);
    subplot(2,2,1);
    semilogy(angle,proc.p11_raw(idx,:))
    xlabel('Angle   (degree)')
    ylabel('P11 Raw Count')
    set(gca,'xlim',[20 145]);
    subplot(2,2,2);
    semilogy(angle,proc.p11_uncal(idx,:))
    xlabel('Angle   (degree)')
    ylabel('P11 Count - Background')
    set(gca,'xlim',[20 145]);
    subplot(2,2,3);
    loglog(proc.Angles,proc.P11(idx,:))
    xlabel('Angle   (degree)')
    ylabel('LISST-VSF P11')
    set(gca,'xlim',[0 160]);
    subplot(2,2,4);
    angle = 15:1:150;
    semilogy(angle,p11_absCal_withbgd,'b.','markersize',10.0)
    hold on
    semilogy(angle,p11_absCal,'r.','markersize',10.0)
    semilogy(angle,P11_lisst(33:168),'k.','markersize',10.0)
    xlabel('Angle   (degree)')
    ylabel('Calibrated P11 (m^-^1)')
    legend('Include background','Substract background','LISST-VSF')
    legend boxoff
    set(gca,'xlim',[20 145]);
 end   
     
 