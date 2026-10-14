function deadline_simulation(mode)
% Run and save the manuscript simulations; default mode uses fixed coefficients.
% Requires Control System Toolbox. Optional redesign also requires Robust Control and Optimization Toolboxes.
if nargin<1,mode='paper';end
if ~any(strcmpi(mode,{'paper','redesign'})),error('Use paper or redesign.');end
run_model_study(strcmpi(mode,'redesign'));
run_ieee_stress_sensitivity();
export_manuscript_tables();
fprintf('Finished. Outputs: %s\n',fullfile(fileparts(mfilename('fullpath')),'results'));
end

function run_model_study(redesign)
% 27 design and 48 evaluation cases; all output trajectories are retained.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
par=struct('L',.419,'ml',.065,'Jl',.0038,'Js',.00208,'bs',.04,'kl',1.3, ...
 'Rm',2.6,'etam',.69,'etag',.90,'kg',70,'km',.00768,'kt',.00768);
[P0,meta0]=ieee_lumped_model(par);po=eig(P0.A);pc=po(imag(po)>1e-5);[~,ix]=max(imag(pc));pflex=pc(ix);
wn_model=abs(pflex);fn_model=wn_model/(2*pi);zeta_model=-real(pflex)/wn_model;wd_model=imag(pflex);fd_model=wd_model/(2*pi);
fn_paper=21.3075;zeta_paper=.3085;wd_paper=fn_paper*sqrt(1-zeta_paper^2);
delay_paper=pi/wd_paper;
comparison=table(fn_model,wn_model,zeta_model,wd_model,fd_model,fn_paper,wd_paper,zeta_paper,delay_paper, ...
 'VariableNames',{'Model_natural_frequency_Hz','Model_natural_omega_rad_s','Model_mode_zeta', ...
 'Model_damped_omega_rad_s','Model_damped_frequency_Hz','IEEE_reported_natural_omega_rad_s', ...
 'IEEE_implied_damped_omega_rad_s','IEEE_reported_zeta','IEEE_shaper_half_period_s'});
writetable(comparison,fullfile(out,'ieee_modal_consistency.csv'));

s=tf('s');G=P0(1,:);wcPID=.70*wd_model;wcHinf=.50*wd_model;
if redesign
 Cpid=pidtune(G,'PIDF',wcPID);
 Ws=(s/2+wcHinf)/(s+.02*wcHinf);Wu=.10;Wt=.5*(s/(.3*wd_model)+.01)/(s/(30*wd_model)+1);
 [Ch,~,gamma]=mixsyn(G,Ws,Wu,Wt);
else
 [Cpid,Ch,fixed,coefmeta]=deadline_fixed_controllers();gamma=coefmeta.hinf_gamma;
end
Tpid=study_close_outputs(P0,Cpid);Th=study_close_outputs(P0,Ch);
assert(isstable(Tpid),'Nominal IEEE plant PIDF loop is unstable');
assert(isstable(Th),'Nominal IEEE plant H-infinity loop is unstable');

Fpaper=struct('a',[.3966 .4296 .1551 .0187],'d',[0 .155 .310 .465],'lp',true);
kappa=exp(real(pflex)*pi/imag(pflex));
Fzvd=struct('a',[1 2*kappa kappa^2]/(1+kappa)^2, ...
 'd',[0 pi/imag(pflex) 2*pi/imag(pflex)],'lp',false);
if ~redesign,Fzvd=fixed.ZVD;Fpaper=fixed.iISC;end
factors=[.75 1 1.25];[fl,fm,fk]=ndgrid(factors,factors,factors);designGrid=[fl(:) fm(:) fk(:)];
ns=size(designGrid,1);
held=[ ...
 1.04568275981162 1.20007521571676 1.0591972546297; ...
 1.19167041280491 0.894122335139818 1.2326527122379; ...
 0.982105778012857 1.18318683687014 1.03225054967906; ...
 1.13147878205541 0.882608780112963 1.04714619715421; ...
 0.794405754548558 1.11383457555753 0.818283675937252; ...
 1.05193784904769 0.831593335969654 0.893238757868595; ...
 0.791169378465167 1.02101567949617 1.00629936930305; ...
 0.837269387515594 0.903436528720904 0.759620486788193; ...
 0.757326004843102 0.867760983190404 0.991493330104688; ...
 1.16295796528131 0.828830301726203 1.02307877232044; ...
 1.06644929821924 0.970440274174861 0.812747035849685; ...
 1.17483569424494 1.12438066720499 0.928376361655957; ...
 0.86877167293452 1.11241786447495 0.945721319555058; ...
 0.756585182999533 1.08361511808362 0.948683061211505; ...
 0.899209293256374 0.951743822026904 0.921973479958472; ...
 1.09421880559442 1.01486393832762 0.772674032379687; ...
 0.953685264458366 1.14236351435028 1.00067264908834; ...
 0.864315530425062 0.830541526447251 0.752714884947609; ...
 0.772791381378108 1.21453908282931 0.938727196382389; ...
 1.00099914296968 1.15672021728589 0.760742973433354; ...
 1.17704672033867 0.797913869438577 1.12739375264545; ...
 0.826039400350022 0.765484866078851 1.17715785463512; ...
 1.01934380504501 1.03204112277943 1.03470890212818; ...
 0.994609363276365 1.05798683339655 1.18802150770759; ...
 0.818227323930354 1.07394475858123 1.12862355032071; ...
 0.839546294293759 1.22579882051896 0.965036884258463; ...
 0.850371191920464 1.07509697506558 0.834502508250038; ...
 1.15761350556037 1.00678899256195 1.04830368432679; ...
 1.07464887329652 1.22808545087435 0.859060712688227; ...
 0.839438851495279 1.02986315153828 0.947366420822732; ...
 1.06857181725047 1.04498675926879 0.807191517168357; ...
 0.795785342495213 1.21923778246239 0.909786395036354; ...
 0.86072001995652 0.980255724251573 1.18263404284974; ...
 0.934816312959206 1.16695676885702 1.17679663126499; ...
 1.0552240692207 0.94632146968611 0.787579190920855; ...
 0.826529658413185 0.924199173849723 0.879378791308633; ...
 1.23857812685521 1.04221009538977 0.881464029187865; ...
 0.918876923663541 1.08645084644012 1.00473044835718; ...
 0.779703070874818 1.0115779763489 0.924219659736757; ...
 0.876040220315324 0.894622936176477 1.08880580696836; ...
 1.03232305336238 1.09796353416117 0.840353593678461; ...
 0.832349342448216 0.884284597112195 0.76508970326041; ...
 0.781653055728789 0.990389132568913 1.24953340204268; ...
 1.16880813395941 1.08410422645021 1.1758648669396; ...
 0.796635233076244 1.17608596841084 1.14686546523962; ...
 1.01009928979015 0.820418111384052 1.23978790943082; ...
 0.932353432865856 1.15672191487718 1.19441344507003; ...
 1.18926252289756 0.94882681299259 0.8286079499388];

pars=[designGrid;held];ncase=size(pars,1);
dt=.001;t=(0:dt:8)';amp=30*pi/180;starts=[0 2.5 5];targets=[amp 1.5*amp 0];Tb=1;
r=zeros(size(t));prev=0;
for j=1:numel(starts)
 r=r+(targets(j)-prev)*study_smoothmove(t,starts(j),Tb);prev=targets(j);
end
delayGrid=0:.025:.5;
if redesign
 nd=numel(delayGrid);designT=cell(1,ns);Ytrain=cell(1,ns);
for j=1:ns
 q=study_variant(par,designGrid(j,:));P=ieee_lumped_model(q);Q=study_close_outputs(P,Cpid);
 designT{j}=Q;Ytrain{j}=lsim(Q,r,t);
end
[aLP,etaLP,flagLP]=study_robust_shaper_lp(Ytrain,t,delayGrid,starts,targets,Tb);

 if flagLP<=0,error('Scenario FIR optimization failed; no paper comparison generated.');end
 FLP=struct('a',aLP,'d',delayGrid,'lp',false);
else
 FLP=fixed.FIR;flagLP=1;etaLP=0.021074654283300843;
end

methods={'PIDF','PIDF with model ZVD','PIDF with IEEE iISC','PIDF with scenario LP','H infinity'};
Cs={Cpid,Cpid,Cpid,Cpid,Ch};Fs={study_identity_filter(),Fzvd,Fpaper,FLP,study_identity_filter()};Ts={Tpid,Tpid,Tpid,Tpid,Th};
if isempty(FLP),methods(4)=[];Cs(4)=[];Fs(4)=[];Ts(4)=[];end
nm=numel(methods);rows={};nominal=cell(1,nm);
nominalCase=find(all(pars==1,2),1);
% responses{case,method}: columns [theta(rad), beta(rad), Vref(V)].
responses=cell(ncase,nm);caseParameters=cell(ncase,1);
deadlineErrors_rad=nan(ncase,nm,3);commands=zeros(numel(t),nm);
for k=1:nm,commands(:,k)=study_shape_reference(t,r,Fs{k});end
for j=1:ncase
 q=study_variant(par,pars(j,:));caseParameters{j}=q;P=ieee_lumped_model(q);
 for k=1:nm
  Tcl=study_close_outputs(P,Cs{k});ok=isstable(Tcl);
  cmd=commands(:,k);
  if ok
   y=lsim(Tcl,cmd,t);[met,deadlineErrors]=study_score(t,r,y,starts,targets,Tb,Fs{k}.d(end),10,dt);
  else
   y=nan(numel(t),3);met=nan(1,9);deadlineErrors=nan(1,3);
  end
  rows(end+1,:)=[{j,study_ternary(j<=ns,'design','heldout'),methods{k},pars(j,1),pars(j,2),pars(j,3),ok},num2cell(met)];
  responses{j,k}=y;
  deadlineErrors_rad(j,k,:)=reshape(deadlineErrors,1,1,3);
  if j==nominalCase,nominal{k}=y;end
 end
 fprintf('IEEE model test case %d/%d complete\n',j,ncase);
end
Tresults=cell2table(rows,'VariableNames',{'Case','Partition','Method','Length_factor','Mass_factor','Stiffness_factor','Stable', ...
 'Angle_NRMSE','Beta_RMS_deg','Beta_peak_deg','Voltage_peak_V','Voltage_rate_peak_V_s','Voltage_limit_pass','Deadline_pass_count','Worst_deadline_error_deg','Saturation_fraction'});
writetable(Tresults,fullfile(out,'ieee_model_comparison.csv'));
writetable(Tresults(Tresults.Case==nominalCase,:),fullfile(out,'Table_nominal_results.csv'));
lpw=zeros(numel(delayGrid),1);if ~isempty(FLP),lpw(1:numel(FLP.a))=FLP.a(:);end
writetable(array2table([delayGrid(:),lpw],'VariableNames',{'Delay_s','Weight'}),fullfile(out,'ieee_lp_shaper.csv'));
Pnom=ieee_lumped_model(par);
meta=struct('par',par,'modelPole',pflex,'wnModel_rad_s',wn_model,'fModel_Hz',fn_model,'wdModel_rad_s',wd_model,'fdModel_Hz',fd_model,'zetaModel',zeta_model, ...
 'reportedOmega',fn_paper,'reportedZeta',zeta_paper,'reportedImpliedWd',wd_paper, ...
 'publishedShaperHalfPeriod',delay_paper,'pid',Cpid,'hinf',Ch,'gamma',gamma, ...
 'delayGrid',delayGrid,'LPweights',FLP,'LPobjective_rad',etaLP,'LPflag',flagLP, ...
 'designCases',designGrid,'heldoutCases',held,'methods',{methods},'dt',dt,'t',t,'reference',r);
meta.nominalCase=nominalCase;meta.caseFactors=pars;
meta.outputNames={'theta','beta','Vref'};meta.outputUnits={'rad','rad','V'};
meta.motionStarts_s=starts;meta.targets_rad=targets;meta.deadlines_s=starts+Tb+.9;
meta.deadlineTolerance_rad=.02*max(abs(targets));meta.matlabVersion=version;
meta.mode='paper';if redesign,meta.mode='redesign';end
save(fullfile(out,'ieee_model_study.mat'),'meta','Tresults','Pnom','Tpid','Th', ...
 'Fs','Cs','Ts','nominal','comparison','responses','commands','caseParameters', ...
 'deadlineErrors_rad','-v7');
fprintf('\nIEEE modal comparison:\n');disp(comparison);
fprintf('\nNominal method summary:\n');disp(Tresults(Tresults.Case==nominalCase,:));
fprintf('\nHeld-out summary:\n');disp(groupsummary(Tresults(strcmp(Tresults.Partition,'heldout'),:),'Method',{'mean','max'}, ...
 {'Angle_NRMSE','Beta_RMS_deg','Beta_peak_deg','Voltage_peak_V','Voltage_limit_pass','Deadline_pass_count','Saturation_fraction'}));
end
function T=study_close_outputs(P,C)
[Ap,Bp,Cp,~]=ssdata(P);[Ak,Bk,Ck,Dk]=ssdata(ss(C));np=size(Ap,1);
cy=Cp(1,:);Ac=[Ap-Bp*Dk*cy Bp*Ck;-Bk*cy Ak];Bc=[Bp*Dk;Bk];
Cc=[Cp(1:2,:) zeros(2,size(Ak,1));-Dk*cy Ck];Dc=[0;0;Dk];
T=ss(Ac,Bc,Cc,Dc);
end
function q=study_variant(p,f)
q=p;q.L=p.L*f(1);q.ml=p.ml*f(2);q.Jl=p.Jl*f(1)^2*f(2);q.kl=p.kl*f(3);
end
function F=study_identity_filter(),F=struct('a',1,'d',0,'lp',false);end
function y=study_shape_reference(t,r,F)
if isempty(F),y=r;return;end
if F.lp,rf=lsim(tf(15,[1 15]),r,t);else,rf=r;end
y=zeros(size(r));for i=1:numel(F.a),y=y+F.a(i)*interp1(t,rf,t-F.d(i),'linear',0);end
end
function v=study_smoothmove(t,st,T)
x=max(0,min(1,(t-st)/T));v=10*x.^3-15*x.^4+6*x.^5;
end
function [a,eta,flag]=study_robust_shaper_lp(Y,t,d,starts,targets,Tb)
% Sampled minimax flexure design with nonnegative FIR taps, voltage envelope, and endpoint tolerance.
nd=numel(d);dt=t(2)-t(1);stride=5;ix=1:stride:numel(t);A=[];b=[];amp=max(abs(targets));
for j=1:numel(Y)
 z=Y{j};B=zeros(numel(t),nd,3);
 for h=1:nd
  sh=round(d(h)/dt);B(sh+1:end,h,:)=z(1:end-sh,:);
 end
 tip=squeeze(B(ix,:,2));volt=squeeze(B(ix,:,3));
 A=[A;tip -ones(size(tip,1),1);-tip -ones(size(tip,1),1);volt zeros(size(volt,1),1);-volt zeros(size(volt,1),1)];
 b=[b;zeros(2*size(tip,1),1);10*ones(2*size(volt,1),1)];
 for e=1:numel(starts)
  td=starts(e)+Tb+.9;[~,it]=min(abs(t-td));row=squeeze(B(it,:,1));target=targets(e);tol=.02*max(abs(targets));
  A=[A;row 0;-row 0];b=[b;target+tol;-target+tol];
 end
end
opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8);
[x,~,flag]=linprog([zeros(nd,1);1],A,b,[ones(1,nd) 0],1,zeros(nd+1,1),[],opts);
if flag>0,a=x(1:nd)';eta=x(end);else,a=[];eta=NaN;end
end
function [met,dead]=study_score(t,r,y,starts,targets,Tb,Tshape,vmax,dt)
% Tracking uses 0-7 s; flexure and voltage metrics use the full 0-8 s record.
% Voltage envelope is screened; no saturation or slew limiter is simulated.
active=t>=0&t<=7;err=y(active,1)-r(active);beta=y(:,2);v=y(:,3);dead=[];
for e=1:numel(starts)
 td=starts(e)+Tb+.9;[~,it]=min(abs(t-td));dead(e)=abs(y(it,1)-targets(e));
end
met=[sqrt(mean(err.^2))/max(abs(r)),sqrt(mean(beta.^2))*180/pi,max(abs(beta))*180/pi, ...
 max(abs(v)),max(abs(diff(v)))/dt,max(abs(v))<=vmax*(1+1e-6),sum(dead<=.02*max(abs(targets))),max(dead)*180/pi,mean(abs(v)>=vmax)];
end
function y=study_ternary(c,a,b),if c,y=a;else,y=b;end,end

function [P,meta]=ieee_lumped_model(par)
% State: [theta, theta_dot, beta, beta_dot]. Input: Vref (V).
% Plant outputs: theta (rad), beta (rad), L*beta (m).
% Model source: Alhassan et al., IEEE Access 11 (2023), 86711-86726.
Js=par.Js; Jl=par.Jl; bs=par.bs; kl=par.kl;
Kv=par.etag*par.kg*par.etam*par.kt/par.Rm;
bm=bs+Kv*par.kg*par.km;
A=[0 1 0 0;
   0 -bm/Js kl/Js 0;
   0 0 0 1;
   0 bm/Js -kl*(1/Js+1/Jl) 0];
B=[0;Kv/Js;0;-Kv/Js];
C=[1 0 0 0;0 0 1 0;0 0 par.L 0];
P=ss(A,B,C,zeros(3,1));
P.InputName={'Vref'}; P.OutputName={'theta','beta','tip_deflection'};
M=[Js+Jl Jl;Jl Jl];
Minv=inv(M);
Aq=[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];
A2=zeros(4); A2(1,2)=1;A2(3,4)=1;
Gq=[0 -bm 0 0;0 0 -kl 0];
Bq=[Kv;0];
A2([2 4],:)=Minv*Gq;
B2=zeros(4,1); B2([2 4])=Minv*Bq;
meta=struct('M',M,'Kv_Nm_per_V',Kv,'equation5_6_A',A,'lagrange_A',A2, ...
 'equation5_6_B',B,'lagrange_B',B2,'A_difference',norm(A-A2,'fro'), ...
 'B_difference',norm(B-B2),'open_loop_poles',eig(A), ...
 'electrical_damping_Nms_per_rad',Kv*par.kg*par.km);
end

function [Cpid, Chinf, shapers, metadata] = deadline_fixed_controllers()
% Full precision manuscript coefficients; controller input is error in rad and output is V.
Cpid = pid(13.826354621552515, 33.880648560318534, 1.0185788286458155, 0.033126913205235989);
Ak = [ ...
    -0.20313685316006516 -1.1250409952694259e-14 -1.350957189983633e-13 -3.3438802354954135e-15 -2.9129578160807874e-14 -2.0638089650070391e-15; ...
    0.0071072561367217912 -609.41087958143305 127.99615621068116 -9.5141217116094217e-05 -0.00082880465959244527 -5.872019420479809e-05; ...
    0.033837788610847339 -0.0015240048287180116 -0.018300357821652503 0.99954703076250528 -0.0039459555602639499 -0.00027956801900119293; ...
    219461.76433874344 -4529.4756393598091 -119831.1662617543 -3042.9388929304387 -25044.040330839431 -1836.1856881786357; ...
    6.1779160398555015e-10 -2.7824436076951829e-11 -3.3411779778124745e-10 -8.2700614690946413e-12 -7.2043071220518415e-11 0.99999999999489575; ...
    -219461.76433873698 4529.4756393594935 119831.16626175081 3042.9388929303523 24701.935067680795 1836.1856881785827];
Bk = [3.9999999999999991; -2.5270754308001559e-05; -0.00012031456673015562; -2.8193413962074575e-11; -2.1966367241502104e-12; 4.7375553508861499e-12];
Ck = [3555.0416161765102 -73.372573331392005 -1941.1344124526254 -48.443179168126115 -415.81050302044605 -29.744208774459242];
Dk = 0;
Chinf = ss(Ak, Bk, Ck, Dk);
shapers.identity = struct('a',1,'d',0,'lp',false);
shapers.ZVD = struct('a',[0.5358022898976017 0.39236576129453571 0.071831948807862714], 'd',[0 0.15465399826344944 0.30930799652689889], 'lp',false);
shapers.FIR = struct('a',[0.36412339161952517 0 0 0 0 0 0 0 0 0 0 0.19031435974162295 0.10996078010927246 0 0 0 0 0 0 0 0.33560146852957945], 'd',[0 0.025000000000000001 0.050000000000000003 0.075000000000000011 0.10000000000000001 0.125 0.15000000000000002 0.17500000000000002 0.20000000000000001 0.22500000000000001 0.25 0.27500000000000002 0.29999999999999999 0.32499999999999996 0.34999999999999998 0.375 0.40000000000000002 0.42499999999999999 0.45000000000000001 0.47499999999999998 0.5], 'lp',false);
shapers.iISC = struct('a',[0.3966 0.4296 0.1551 0.0187], 'd',[0 0.155 0.310 0.465], 'lp',true);
metadata.wn_rad_s = 21.327236479382776;
metadata.wd_rad_s = 20.313685316031496;
metadata.zeta = 0.30461299110912615;
metadata.hinf_gamma = 1.3594991786584194;
metadata.MATLAB_version = '25.2.0.3123386 (R2025b) Update 3';
end

function run_ieee_stress_sensitivity()
% Six scenarios with fixed coefficients; controller initial states are zero.
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
D=load(fullfile(out,'ieee_model_study.mat'),'meta','Pnom','Cs','Fs');
meta=D.meta; par=meta.par; t=meta.t(:); r=meta.reference(:); dt=meta.dt;
starts=[0 2.5 5]; targets=[30 45 0]*pi/180; Tb=1;
methods=meta.methods(:); rows={};
scenarios=struct('name',{},'Lfactor',{},'massFactor',{},'stiffnessFactor',{},'dampingFactor',{},'beta0_deg',{});
scenarios(1)=struct('name','Nominal parameters','Lfactor',1,'massFactor',1,'stiffnessFactor',1,'dampingFactor',1,'beta0_deg',0);
scenarios(2)=struct('name','Soft link, half stiffness','Lfactor',1,'massFactor',1,'stiffnessFactor',.5,'dampingFactor',1,'beta0_deg',0);
scenarios(3)=struct('name','Reduced hub damping','Lfactor',1,'massFactor',1,'stiffnessFactor',1,'dampingFactor',.5,'beta0_deg',0);
scenarios(4)=struct('name','Increased link mass','Lfactor',1,'massFactor',1.5,'stiffnessFactor',1,'dampingFactor',1,'beta0_deg',0);
scenarios(5)=struct('name','Initial flexure 2 deg','Lfactor',1,'massFactor',1,'stiffnessFactor',1,'dampingFactor',1,'beta0_deg',2);
scenarios(6)=struct('name','Combined soft, heavy, low damping','Lfactor',1,'massFactor',1.25,'stiffnessFactor',.5,'dampingFactor',.5,'beta0_deg',0);
responses=cell(numel(scenarios),numel(methods));caseParameters=cell(numel(scenarios),1);
deadlineErrors_deg=nan(numel(scenarios),numel(methods),3);
commands=zeros(numel(t),numel(methods));
for k=1:numel(methods),commands(:,k)=stress_shape_reference(t,r,D.Fs{k});end
for j=1:numel(scenarios)
    sc=scenarios(j); q=par;
    q.L=par.L*sc.Lfactor; q.ml=par.ml*sc.massFactor;
    q.Jl=par.Jl*sc.Lfactor^2*sc.massFactor;
    q.kl=par.kl*sc.stiffnessFactor; q.bs=par.bs*sc.dampingFactor;
    caseParameters{j}=q;
    P=ieee_lumped_model(q);
    for k=1:numel(methods)
        Tcl=stress_close_outputs(P,D.Cs{k});
        stable=isstable(Tcl);
        cmd=commands(:,k);
        if stable
            [A,~,~,~]=ssdata(Tcl); x0=zeros(size(A,1),1);
            x0(3)=sc.beta0_deg*pi/180; % Initial relative flexure angle.
            y=lsim(Tcl,cmd,t,x0);
            active=t<=7; err=y(active,1)-r(active);
            deadline=zeros(1,3);
            for e=1:3
                td=starts(e)+Tb+.9; [~,it]=min(abs(t-td));
                deadline(e)=abs(y(it,1)-targets(e))*180/pi;
            end
            met=[sqrt(mean(err.^2))/max(abs(r)),sqrt(mean(y(:,2).^2))*180/pi, ...
                max(abs(y(:,2)))*180/pi,max(abs(y(:,3))), ...
                sum(deadline<=.9),max(deadline)];
        else
            met=nan(1,6);y=nan(numel(t),3);deadline=nan(1,3);
        end
        responses{j,k}=y;
        deadlineErrors_deg(j,k,:)=reshape(deadline,1,1,3);
        rows(end+1,:)=[{sc.name,sc.Lfactor,sc.massFactor,sc.stiffnessFactor, ...
            sc.dampingFactor,sc.beta0_deg,methods{k},stable},num2cell(met)];
    end
end
T=cell2table(rows,'VariableNames',{'Scenario','Length_factor','Mass_factor', ...
 'Stiffness_factor','Hub_damping_factor','Initial_beta_deg','Method','Stable', ...
 'Angle_NRMSE','Beta_RMS_deg','Beta_peak_deg','Voltage_peak_V', ...
 'Deadline_passes','Worst_deadline_error_deg'});
writetable(T,fullfile(out,'ieee_stress_sensitivity.csv'));
Cs=D.Cs;Fs=D.Fs;
save(fullfile(out,'ieee_stress_sensitivity.mat'),'T','scenarios','meta','Cs','Fs', ...
 'responses','commands','caseParameters','deadlineErrors_deg','-v7');
disp(T);
end
function T=stress_close_outputs(P,C)
[Ap,Bp,Cp,~]=ssdata(P); [Ak,Bk,Ck,Dk]=ssdata(ss(C));
cy=Cp(1,:); Ac=[Ap-Bp*Dk*cy Bp*Ck;-Bk*cy Ak]; Bc=[Bp*Dk;Bk];
Cc=[Cp(1:2,:) zeros(2,size(Ak,1));-Dk*cy Ck]; Dc=[0;0;Dk];
T=ss(Ac,Bc,Cc,Dc);
end
function y=stress_shape_reference(t,r,F)
if F.lp, rf=lsim(tf(15,[1 15]),r,t); else, rf=r; end
y=zeros(size(r));
for i=1:numel(F.a), y=y+F.a(i)*interp1(t,rf,t-F.d(i),'linear',0); end
end

function export_manuscript_tables()
% Export evaluation and sensitivity summaries from the simulated metrics.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results');
R=readtable(fullfile(out,'ieee_model_comparison.csv'));
H=R(strcmp(R.Partition,'heldout'),:);
methods=unique(string(H.Method),'stable');rows=cell(numel(methods),15);
for k=1:numel(methods)
    q=H(strcmp(string(H.Method),methods(k)),:);
    rows(k,:)={char(methods(k)),height(q),sum(q.Stable),height(q), ...
        mean(q.Angle_NRMSE,'omitnan'),std(q.Angle_NRMSE,'omitnan'), ...
        mean(q.Beta_RMS_deg,'omitnan'),std(q.Beta_RMS_deg,'omitnan'), ...
        mean(q.Beta_peak_deg,'omitnan'),std(q.Beta_peak_deg,'omitnan'), ...
        max(q.Beta_peak_deg,[],'omitnan'),max(q.Voltage_peak_V,[],'omitnan'), ...
        max(q.Voltage_rate_peak_V_s,[],'omitnan'),sum(q.Deadline_pass_count,'omitnan'),3*height(q)};
end
T=cell2table(rows,'VariableNames',{'Method','Heldout_cases','Stable_cases','Total_cases', ...
    'Mean_NRMSE','SD_NRMSE','Mean_beta_RMS_deg','SD_beta_RMS_deg', ...
    'Mean_beta_peak_deg','SD_beta_peak_deg','Worst_beta_peak_deg','Worst_voltage_V', ...
    'Worst_voltage_slew_V_s','Deadline_passes','Deadline_checks'});
writetable(T,fullfile(out,'ieee_heldout_summary.csv'));
writetable(T(:,[1:3 5:15]),fullfile(out,'Table_heldout_results.csv'));
S=readtable(fullfile(out,'ieee_stress_sensitivity.csv'));
rows=cell(numel(methods),6);
for k=1:numel(methods)
    q=S(strcmp(string(S.Method),methods(k)),:);
    rows(k,:)={char(methods(k)),height(q),sum(q.Stable), ...
        sum(q.Deadline_passes,'omitnan'),3*height(q),max(q.Voltage_peak_V,[],'omitnan')};
end
Q=cell2table(rows,'VariableNames',{'Method','Scenarios','Stable_scenarios', ...
    'Deadline_passes','Deadline_checks','Worst_voltage_V'});
writetable(Q,fullfile(out,'stress_deadline_summary.csv'));
stressMethods=unique(string(S.Method));rows=cell(numel(stressMethods),11);
for k=1:numel(stressMethods)
    q=S(strcmp(string(S.Method),stressMethods(k)),:);
    rows(k,:)={char(stressMethods(k)),height(q),sum(q.Stable), ...
        mean(q.Angle_NRMSE,'omitnan'),max(q.Angle_NRMSE,[],'omitnan'), ...
        mean(q.Beta_RMS_deg,'omitnan'),max(q.Beta_peak_deg,[],'omitnan'), ...
        max(q.Voltage_peak_V,[],'omitnan'),sum(q.Deadline_passes,'omitnan'), ...
        3*height(q),max(q.Worst_deadline_error_deg,[],'omitnan')};
end
Q=cell2table(rows,'VariableNames',{'Method','Scenarios','Stable_scenarios', ...
    'Mean_NRMSE','Max_NRMSE','Mean_beta_RMS_deg','Worst_beta_peak_deg', ...
    'Worst_voltage_V','Deadline_passes','Deadline_checks','Worst_deadline_error_deg'});
writetable(Q,fullfile(out,'Table_ieee_stress_summary.csv'));
end
