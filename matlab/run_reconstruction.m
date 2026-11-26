%% QP-Based Vehicle Speed Planning -- V1.0 reconstruction
% Reconstructed from the October 2026 report, equations (7)-(27).
% Requires MATLAB Optimization Toolbox; R2020a or later for exportgraphics.
% All model parameters are fixed to the documented baseline.
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'results','matlab');
if ~exist(out,'dir'), mkdir(out); end

%% 1. Spatial grid and decision-vector indices
L = 60; ds = 0.5; s = (0:ds:L)';
N = numel(s); M = N-1; n = N+M;
stopS = [0 20 40 60];
stopIdx = round(stopS/ds)+1; % MATLAB uses one-based indexing.
ib = 1:N; ia = N+(1:M);
vRef = 5; vMax = 6; aMax = 1.5; aMin = -2;
wTrack = 0.08; weights = [0 1 10 100]; activeTol = 1e-4;
betaRef = vRef^2*ones(N,1);

%% 2. Dynamics, stop equalities and physical bounds
Aeq = spalloc(M+numel(stopIdx),n,3*M+numel(stopIdx));
for i = 1:M
    Aeq(i,ib(i)) = -1;
    Aeq(i,ib(i+1)) = 1;
    Aeq(i,ia(i)) = -2*ds;
end
for k = 1:numel(stopIdx)
    Aeq(M+k,ib(stopIdx(k))) = 1;
end
beq = zeros(size(Aeq,1),1);
lb = [zeros(N,1); aMin*ones(M,1)];
ub = [vMax^2*ones(N,1); aMax*ones(M,1)];
assert(n==241 && isequal(size(Aeq),[124 241]));
assert(rank(full(Aeq))==124 && nnz(Aeq)==364);

%% 3. First difference: D*alpha = alpha(i+1)-alpha(i)
D = spalloc(M-1,M,2*(M-1));
for i = 1:M-1
    D(i,i) = -1; D(i,i+1) = 1;
end
assert(isequal(size(D),[119 120]));
f = [-2*wTrack*betaRef; zeros(M,1)];
constant = wTrack*(betaRef'*betaRef); % Omitted by quadprog.
opts = optimoptions('quadprog','Algorithm','interior-point-convex', ...
    'Display','off','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);

%% 4. Solve each documented weight with the original solver
K = numel(weights); V = zeros(N,K); A = zeros(M,K); activity = zeros(K,M);
summary = zeros(K,19); runData = struct([]);
for k = 1:K
    ws = weights(k);
    H = blkdiag(2*wTrack*speye(N),2*ws*(D'*D));
    assert(isequal(size(H),[241 241]) && norm(H-H','inf')<1e-12);
    assert(min(eig(full(H))) >= -1e-9); % PSD, not necessarily PD.
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],opts);
    assert(exitflag>0 && numel(x)==n,'quadprog failed at ws=%g',ws);
    beta = x(ib); alpha = x(ia);
    % Clamp only for physical-speed display, never for KKT or costs.
    velocity = sqrt(max(beta,0));
    V(:,k)=velocity; A(:,k)=alpha;
    rawTrack = sum((beta-betaRef).^2);
    rawSmooth = sum((D*alpha).^2);
    weightedTrack = wTrack*rawTrack; weightedSmooth = ws*rawSmooth;
    J = weightedTrack+weightedSmooth;
    objResidual = abs(J-(fval+constant));
    assert(objResidual<1e-7*max(1,J),'Objective reconstruction mismatch');
    % 5. Check all four KKT conditions using the UNMODIFIED solver output.
    eqInf = norm(Aeq*x-beq,inf);
    boundViolation = max([0;lb-x;x-ub]);
    statInf = norm(H*x+f+Aeq'*lambda.eqlin-lambda.lower+lambda.upper,inf);
    dualViolation = max([0;-lambda.lower;-lambda.upper]);
    compInf = max(norm(lambda.lower.*(x-lb),inf),norm(lambda.upper.*(ub-x),inf));
    kktPass = max([eqInf,boundViolation,statInf,dualViolation,compInf]) <= 1e-7;
    activeAcc = abs(alpha-aMax)<activeTol;
    activeBrake = abs(alpha-aMin)<activeTol;
    activity(k,:) = (double(activeAcc)-double(activeBrake))';
    stopInf = norm(beta(stopIdx),inf);
    summary(k,:) = [ws rawTrack weightedTrack rawSmooth weightedSmooth J ...
        sum(activeAcc) sum(activeBrake) sum(activeAcc|activeBrake) ...
        eqInf boundViolation statInf dualViolation compInf exitflag ...
        kktPass fval objResidual stopInf];
    prefix = sprintf('ws_%g_',ws);
    writetable(table((0:N-1)',s,beta,velocity,'VariableNames', ...
        {'node','s_m','beta','speed_mps'}),fullfile(out,[prefix 'nodes.csv']));
    writetable(table((0:M-1)',s(1:end-1),s(2:end),alpha,activeAcc,activeBrake, ...
        'VariableNames',{'segment','s_start_m','s_end_m','acceleration_mps2', ...
        'active_upper','active_lower'}),fullfile(out,[prefix 'segments.csv']));
    writetable(table((0:n-1)',lambda.lower,lambda.upper,'VariableNames', ...
        {'variable_index','lambda_lower','lambda_upper'}),fullfile(out,[prefix 'multipliers.csv']));
    writetable(table((0:size(Aeq,1)-1)',lambda.eqlin,'VariableNames', ...
        {'row_index','nu'}),fullfile(out,[prefix 'equality_multipliers.csv']));
    runData(k).ws=ws; runData(k).x=x; runData(k).H=H;
    runData(k).lambda=lambda; runData(k).output=output;
    fprintf('ws=%g | tracking %.8f | smooth %.8f | active %d | KKT %.3e\n', ...
        ws,weightedTrack,rawSmooth,sum(activeAcc|activeBrake), ...
        max([eqInf,boundViolation,statInf,dualViolation,compInf]));
    if ~kktPass, warning('KKT audit threshold exceeded for ws=%g',ws); end
end
names = {'ws','raw_tracking','weighted_tracking','raw_smoothness','weighted_smoothness', ...
    'J','active_acceleration','active_braking','active_total','equality_inf', ...
    'bound_violation','stationarity_inf','dual_violation','complementarity_inf', ...
    'exitflag','kkt_pass','quadprog_fval','objective_residual','stop_beta_inf'};
T = array2table(summary,'VariableNames',names);
writetable(T,fullfile(out,'summary.csv')); disp(T);
save(fullfile(out,'workspace.mat'),'runData','Aeq','beq','D','f','lb','ub','s','weights','T');
fid=fopen(fullfile(out,'environment.txt'),'w');
fprintf(fid,'MATLAB: %s\nPlatform: %s\nAeq rank: %d\n',version,computer,rank(full(Aeq)));
fclose(fid);

%% 6. Plot the documented experiments (no additional optimization features)
labels = compose('w_s = %g',weights);
fig=figure('Color','w','Position',[100 100 900 600]);
tiledlayout(2,1);
nexttile; hold on; grid on;
for k=1:K, plot(s,V(:,k),'LineWidth',1.5); end
legend(labels,'Location','best'); yline(vRef,':','HandleVisibility','off');
for st=stopS, xline(st,':','HandleVisibility','off'); end
xlim([0 L]); ylabel('Speed (m/s)'); title('V1.0 reconstruction: MATLAB quadprog');
nexttile; hold on; grid on;
for k=1:K, stairs(s,[A(:,k); A(end,k)],'LineWidth',1.5); end
for st=stopS, xline(st,':','HandleVisibility','off'); end
yline(aMax,':','HandleVisibility','off'); yline(aMin,':','HandleVisibility','off');
xlim([0 L]); xlabel('Path distance (m)'); ylabel('Acceleration (m/s^2)');
exportgraphics(fig,fullfile(out,'profiles.png'),'Resolution',180);
fig=figure('Color','w','Position',[100 100 1000 400]); tiledlayout(1,2);
ax=nexttile; imagesc(ax,s(1:end-1)+ds/2,1:K,activity,[-1 1]);
colormap(ax,[.84 .36 .31;.90 .92 .94;.09 .56 .49]); cb=colorbar(ax,'Ticks',[-1 0 1],'TickLabels',{'Braking','Inactive','Acceleration'});
yticks(1:K); yticklabels(string(weights)); xlabel('Distance (m)'); ylabel('Smoothing weight'); title('Active bounds by segment');
nexttile; bar(1:K,[T.active_acceleration T.active_braking]); xticks(1:K); xticklabels(string(weights));
legend('Acceleration','Braking'); ylabel('Active segments'); xlabel('Smoothing weight'); title('Active counts');
exportgraphics(fig,fullfile(out,'active_set.png'),'Resolution',180);
fig=figure('Color','w','Position',[100 100 900 400]); tiledlayout(1,2);
nexttile; bar(1:K,[T.weighted_tracking T.weighted_smoothness],'stacked');
legend('Weighted tracking','Weighted smoothing'); xticks(1:K); xticklabels(string(weights));
xlabel('Smoothing weight'); ylabel('Cost'); title('Weighted objective');
nexttile; bar(1:K,T.raw_smoothness); xticks(1:K); xticklabels(string(weights));
xlabel('Smoothing weight'); ylabel('Raw smoothness'); title('Acceleration differences');
exportgraphics(fig,fullfile(out,'objective.png'),'Resolution',180);

%% 7. Cross-check independent reconstruction without tuning the model
baselinePath=fullfile(root,'results','independent','summary.csv');
if isfile(baselinePath)
    B=readtable(baselinePath); assert(isequal(T.ws,B.ws));
    C=table(T.ws,T.raw_tracking-B.raw_tracking,T.weighted_tracking-B.weighted_tracking, ...
        T.raw_smoothness-B.raw_smoothness,T.active_total-B.active_total, ...
        'VariableNames',{'ws','delta_raw_tracking','delta_weighted_tracking','delta_raw_smoothness','delta_active_total'});
    writetable(C,fullfile(out,'comparison_to_independent.csv'));
    disp(C);
    % Metric agreement tolerance chosen for the reconstruction audit, not a model parameter.
    assert(all(abs(C.delta_weighted_tracking)<1e-3 & abs(C.delta_raw_smoothness)<1e-3 & C.delta_active_total==0), ...
        'MATLAB differs from independent reference. Investigate solver and active tolerance.');
end
assert(all(T.kkt_pass==1),'At least one MATLAB KKT audit did not pass.');
