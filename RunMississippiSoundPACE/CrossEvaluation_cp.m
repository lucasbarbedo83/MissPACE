function [cp517] = CrossEvaluation_cp(PACE,inwater)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here


fig1=figure;
ax1=axes('box','on')%,'XScale','log','YScale','log')
hold on
plot([0 35],[0 35],'LineStyle','-','Color','k','LineWidth',0.5)
for nst= 1: length(inwater)
  if ~isempty(ACS(nst))&& ~isempty(LISSTVSPACEF(nst))
        plot( median(ACS(nst).cp,2,"omitnan"),median(PACE(nst).proc.cp,2,"omitnan"),...
            'o', 'LineWidth',0.8,'color','k')
        plot( [mean(ACS(nst).cp,2,"omitnan")-std(ACS(nst).cp) ...
            mean(ACS(nst).cp,2,"omitnan")+std(ACS(nst).cp)],...
            [mean(LISSTVSF(nst).cp,2,"omitnan")  mean(LISSTVSF(nst).cp,2,"omitnan")],...
            '-', 'LineWidth',0.8,'color','k')
        plot( [mean(ACS(nst).cp,2,"omitnan")  mean(ACS(nst).cp,2,"omitnan") ],...
            [mean(LISSTVSF(nst).cp,2,"omitnan")-std(LISSTVSF(nst).cp) ...
            mean(LISSTVSF(nst).cp,2,"omitnan")+std(LISSTVSF(nst).cp)],...
            '-', 'LineWidth',0.8,'color','k')
        str=[nstation(nst).name];
        L=str=='_';
        str(L)=' ';

        text(  mean(ACS(nst).cp,2,"omitnan") ,...
            mean(LISSTVSF(nst).cp,2,"omitnan")-std(LISSTVSF(nst).cp) ,...
            str,'fontsize',7,'Rotation',-30,...
            'color','k');
  end
end
axis equal;

set(ax1,'box','on','XLim',[5 20],'YLim',[5 20])
fontSize=12;
grid minor
fontSize=12;
xlabel('ac-s: $c_p(517)$ ($\mu\pm\sigma$ in m$^{-1}$)',...
    'FontSize',fontSize,'FontWeight','Bold','Interpreter','latex')
ylabel('LISST-VSF: $c_p(517)$ ($\mu\pm\sigma$ in m$^{-1}$)',...
    'FontSize',fontSize,'FontWeight','Bold','Interpreter','latex')

hold off
str=['scattering_LISSTVSFvsACS_cp.png'];
print(fig1, '-dpng',str)



end
end