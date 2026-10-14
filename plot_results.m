function plot_results()
% Read saved results and export Figs. 4-7 without rerunning simulation.
root=fileparts(mfilename('fullpath'));
if ~isfile(fullfile(root,'results','ieee_model_study.mat'))
    error('Run deadline_simulation before plot_results.');
end
refresh_ieee_nominal_figure_readable();
export_ieee_heldout_scatter_readable();
refresh_ieee_figures_6_7_readable();
fprintf('Figures saved in %s\n',fullfile(root,'results'));
end

function refresh_ieee_nominal_figure_readable()
root=fileparts(mfilename('fullpath'));
R=load(fullfile(root,'results','ieee_model_study.mat'),'meta','nominal');
t=R.meta.t; r=R.meta.reference; methods=R.meta.methods;
f=figure('Visible','off','Color','w','Position',[80 80 1400 1000]);
tl=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
ax=gobjects(3,1);
for j=1:3
    ax(j)=nexttile(tl,j); hold(ax(j),'on'); grid(ax(j),'on');
    set(ax(j),'FontSize',16,'LineWidth',1.1,'FontName','Arial');
end
for k=1:numel(methods)
    yc=R.nominal{k};
    plot(ax(1),t,yc(:,1)*180/pi,'LineWidth',2.2,'DisplayName',methods{k});
    plot(ax(2),t,yc(:,2)*180/pi,'LineWidth',2.2,'DisplayName',methods{k});
    plot(ax(3),t,yc(:,3),'LineWidth',2.2,'DisplayName',methods{k});
end
plot(ax(1),t,r*180/pi,'k--','LineWidth',2.0,'DisplayName','Reference');
ylabel(ax(1),'Hub angle (deg)','FontSize',19);
ylabel(ax(2),'Deflection angle \beta (deg)','FontSize',19);
yline(ax(3),10,'k:','LineWidth',1.4,'HandleVisibility','off');
yline(ax(3),-10,'k:','LineWidth',1.4,'HandleVisibility','off');
ylabel(ax(3),{'Motor voltage command','V_{ref} (V)'},'FontSize',19);
xlabel(ax(3),'Time (s)','FontSize',19);
for j=1:3
    lg=legend(ax(j),'Location','eastoutside','FontSize',14);
    lg.ItemTokenSize=[36 18];
end
sgtitle(tl,'IEEE two coordinate model with published parameters','FontSize',19);
base=fullfile(root,'results','Fig4_nominal_responses');
exportgraphics(f,[base '.png'],'Resolution',300);
exportgraphics(f,[base '.pdf'],'ContentType','vector');
savefig(f,[base '.fig']);
close(f);
end

function export_ieee_heldout_scatter_readable()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
T=readtable(fullfile(out,'ieee_model_comparison.csv'));
T=T(strcmp(T.Partition,'heldout') & T.Stable==1,:);
methods={'PIDF','PIDF with model ZVD','PIDF with IEEE iISC','PIDF with scenario LP','H infinity'};
labels={'PIDF','PIDF + ZVD','PIDF + IEEE iISC','PIDF + scenario FIR','H infinity'};
f=figure('Visible','off','Color','w','Position',[100 100 1200 750]);
hold on; colors=lines(numel(methods));
for k=1:numel(methods)
    q=T(strcmp(T.Method,methods{k}),:);
    scatter(q.Angle_NRMSE,q.Beta_peak_deg,34,'filled','MarkerFaceColor',colors(k,:), ...
        'MarkerEdgeColor',colors(k,:),'DisplayName',labels{k});
end
set(gca,'FontSize',16,'FontName','Arial','LineWidth',1.1);
grid on;
xlabel('Normalized angle tracking RMSE','FontSize',19);
ylabel('Peak flexure angle (deg)','FontSize',19);
title('Held out parameter cases using the IEEE model','FontSize',19);
legend('Location','eastoutside','Interpreter','none','FontSize',14);
base=fullfile(out,'Fig5_evaluation_cases');
exportgraphics(f,[base '.png'],'Resolution',300);
exportgraphics(f,[base '.pdf'],'ContentType','vector');
savefig(f,[base '.fig']); close(f);
end

function refresh_ieee_figures_6_7_readable()
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results');
names={'PIDF','PIDF with model ZVD','PIDF with IEEE iISC','PIDF with scenario LP','H infinity'};
labels={'PIDF','PIDF + ZVD','PIDF + IEEE iISC','PIDF + scenario FIR','H infinity'};
T=readtable(fullfile(out,'Table_heldout_results.csv'));
[~,ix]=ismember(names,T.Method);T=T(ix,:);
f=figure('Visible','off','Color','w','Position',[100 100 1300 800]);hold on;
colors=lines(numel(names));
for k=1:numel(names)
    errorbar(T.Mean_NRMSE(k),T.Mean_beta_peak_deg(k),T.SD_beta_peak_deg(k), ...
        T.SD_beta_peak_deg(k),T.SD_NRMSE(k),T.SD_NRMSE(k),'o', ...
        'Color',colors(k,:),'MarkerFaceColor',colors(k,:),'LineWidth',1.5,'MarkerSize',8);
    dy=0;if k==2,dy=.065;elseif k==5,dy=-.065;end
    short=labels{k};
    text(T.Mean_NRMSE(k)+.006,T.Mean_beta_peak_deg(k)+dy,short, ...
        'FontSize',14,'Color',colors(k,:));
end
fig67_style_axes(gca);grid on;xlim([.015 .22]);ylim([.1 1.65]);
xlabel('Mean normalized angle tracking RMSE','FontSize',19);
ylabel('Mean peak flexure angle (deg)','FontSize',19);
title('Tracking and vibration trade off across 48 held out cases','FontSize',19);
subtitle('Error bars show one standard deviation across cases','FontSize',16);
fig67_export_figure(f,fullfile(out,'Fig6_mean_tradeoff'));

T=readtable(fullfile(out,'ieee_stress_sensitivity.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 1800 900]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax1=nexttile(tl);hold(ax1,'on');
for k=1:numel(names)
    q=T(strcmp(T.Method,names{k}) & strcmp(T.Scenario,'Combined soft, heavy, low damping'),:);
    scatter(ax1,q.Angle_NRMSE,q.Beta_peak_deg,75,'filled','DisplayName',names{k});
end
fig67_style_axes(ax1);grid(ax1,'on');
xlabel(ax1,'Normalized angle tracking RMSE','FontSize',19);
ylabel(ax1,'Peak flexure angle (deg)','FontSize',19);
title(ax1,'Combined numerical stress case','FontSize',19);
legend(ax1,'Location','best','Interpreter','none','FontSize',14);
ax2=nexttile(tl);hold(ax2,'on');
q=T(strcmp(T.Scenario,'Initial flexure 2 deg'),:);
[~,ix]=ismember(names,q.Method);q=q(ix,:);
bar(ax2,categorical(names,names),q.Beta_RMS_deg);
fig67_style_axes(ax2);grid(ax2,'on');xtickangle(ax2,30);
ylabel(ax2,'RMS flexure angle (deg)','FontSize',19);
title(ax2,{'Whole record RMS with','2 degree initial deflection'},'FontSize',19);
sgtitle(tl,'IEEE model sensitivity tests with fixed controller and shaper coefficients','FontSize',19);
fig67_export_figure(f,fullfile(out,'Fig7_sensitivity'));
end

function fig67_style_axes(ax)
set(ax,'FontSize',16,'FontName','Arial','LineWidth',1.1);
end

function fig67_export_figure(f,base)
exportgraphics(f,[base '.png'],'Resolution',300);
exportgraphics(f,[base '.pdf'],'ContentType','vector');
savefig(f,[base '.fig']);close(f);
end
