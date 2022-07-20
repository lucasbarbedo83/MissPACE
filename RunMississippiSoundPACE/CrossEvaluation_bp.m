function [bp] = CrossEvaluation_bp(PACE,inwater)
cont=0;
test=struct2cell(inwater);

[li,co,le]=size(test);
ccolor=parula(le);
for nst=1:le
    test1=sum(~isnan(inwater(nst).ap(:)));
    test2=sum(~isnan(PACE(nst).proc.p11(:)));
    if test1>10 && test2>10
        cont=cont+1;
        cp517=interp1(inwater(nst).wl_c,median(inwater(nst).cp,1,'omitnan'),517);
        ap517=interp1(inwater(nst).wl_a,median(inwater(nst).ap,1,'omitnan'),517);
        bp517(cont)=cp517-ap517;

        theta=[1:180];
        Tc= mean(inwater(nst).T);
        SA=mean(inwater(nst).S);
        [betasw,bsw]= betasw_HZP2019(517,theta,Tc,SA);

        beambp(cont)=median(PACE(nst).proc.beamb,2,'omitnan')-bsw;
        cramp(cont,:)=ccolor(nst,:);
    end

end

fig1=figure;
ax1=axes('box','on');%,'XScale','log','YScale','log')
hold on
ax1.ColorOrder = cramp;
for i=1:length(cramp)
 plot(ax1,bp517(i),beambp(i),'o', 'LineWidth',1.8,'Color',cramp(i,:))
end
plot([0 35],[0 35],'LineStyle','-','LineWidth',0.5,'color','k')
axis equal;
set(ax1,'box','on','XLim',[0 20],'YLim',[0 20])
fontSize=12;
grid minor
fontSize=12;
xlabel('ac-s: $c_p(517)$ ($\mu\pm\sigma$ in m$^{-1}$)',...
    'FontSize',fontSize,'FontWeight','Bold','Interpreter','latex')
ylabel('LISST-VSF: $c_p(517)$ ($\mu\pm\sigma$ in m$^{-1}$)',...
    'FontSize',fontSize,'FontWeight','Bold','Interpreter','latex')

hold off
%str=['scattering_LISSTVSFvsACS_cp.png'];
%print(fig1, '-dpng',str)

bp=[beambp; bp517]';
end







        %         plot( [mean(ACS(nst).cp,2,"omitnan")-std(ACS(nst).cp) ...
        %             mean(ACS(nst).cp,2,"omitnan")+std(ACS(nst).cp)],...
        %             [mean(LISSTVSF(nst).cp,2,"omitnan")  mean(LISSTVSF(nst).cp,2,"omitnan")],...
        %             '-', 'LineWidth',0.8,'color','k')
        %         plot( [mean(ACS(nst).cp,2,"omitnan")  mean(ACS(nst).cp,2,"omitnan") ],...
        %             [mean(LISSTVSF(nst).cp,2,"omitnan")-std(LISSTVSF(nst).cp) ...
        %             mean(LISSTVSF(nst).cp,2,"omitnan")+std(LISSTVSF(nst).cp)],...
        %             '-', 'LineWidth',0.8,'color','k')
        % str=[nstation(nst).name];
        % L=str=='_';
        %str(L)=' ';

        %text(  mean(ACS(nst).cp,2,"omitnan") ,...
        %   mean(LISSTVSF(nst).cp,2,"omitnan")-std(LISSTVSF(nst).cp) ,...
        %  str,'fontsize',7,'Rotation',-30,...
        %  'color','k');