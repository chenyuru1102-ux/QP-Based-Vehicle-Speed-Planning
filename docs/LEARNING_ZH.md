# 繁體中文學習手冊：從數學推導到 MATLAB

目標不是背程式，而是每段都能自己寫出來、解釋為什麼這樣建模，並能辨認結果是否合理。按下列順序練習，先不要加 V2.0 功能。

## 第 1 課：為什麼不用速度 v 直接當變數？

從 a=dv/dt 出發，鏈鎖律給出 a=(dv/ds)(ds/dt)=v(dv/ds)。又 d(v²)/ds=2v(dv/ds)=2a。定義 beta=v²，因此 d(beta)/ds=2a。假設每個空間區段內加速度 alpha 為常數，積分得到 beta(i+1)-beta(i)=2*ds*alpha(i)。在 ds=0.5 時，係數恰好為 1；不能把這個簡化套到其他 ds。

手算：beta=9、alpha=1、ds=0.5，下一個 beta=10，速度 sqrt(10)=3.1623 m/s。停車點附近的空間式可由積分關係描述，不能直接在 v=0 時除以 v。

## 第 2 課：決策向量與等式限制

60/0.5=120 個區段，但有 121 個節點（包括兩端）。121 個 beta 加 120 個 alpha 得 241 個未知數。每段一條運動等式，共 120 條；另有 4 個停車等式，合計 124 條。MATLAB 從 1 開始編號，停止點的 beta 欄是 [1,41,81,121]，而非 [0,40,80,120]。

```matlab
s = (0:0.5:60)';
N = numel(s); M = N-1;
ib = 1:N; ia = N+(1:M);
Aeq = sparse(M+4,N+M);
for i = 1:M
    Aeq(i,ib(i)) = -1;
    Aeq(i,ib(i+1)) = 1;
    Aeq(i,ia(i)) = -2*0.5;
end
for k = 1:4
    idx = [1 41 81 121];
    Aeq(M+k,idx(k)) = 1;
end
size(Aeq)
nnz(Aeq)
rank(full(Aeq))
```

應分別得到 [124 241]、364、124。第一列非零位置是第 1、2、122 欄；第 121 列只在第 1 欄為 1。右側 beq=zeros(124,1)。界限 lb=[zeros(N,1); -2*ones(M,1)]，ub=[36*ones(N,1); 1.5*ones(M,1)]。

## 第 3 課：目標函數如何變成 H 和 f

E_track=||beta-betaRef||²，betaRef=25*ones(121,1)。展開為 beta' beta - 2 betaRef' beta + betaRef' betaRef。最後一項是常數，不影響最小點。

E_smooth=||D alpha||²=alpha' D' D alpha。D 每列只有 -1、+1，取相鄰加速度差。先練習 4 個 alpha，例如 [1;1;0;-2]，D alpha=[0;-1;-2]，平滑成本=5。

quadprog 的標準式為 0.5*x'*H*x+f'*x。故 H=2*blkdiag(wTrack*I,ws*D'*D)，f=[-2*wTrack*betaRef;zeros(M,1)]。乘 2 是因為標準式前面有 0.5。常數是 0.08*121*625=6050；完整 J 等於 quadprog fval+6050。不可直接把 fval 當完整成本。

H 為半正定：對任意向量 y，y'D'Dy=||Dy||² >= 0。ws=0 時 alpha block 全零，但 alpha 由 beta 唯一決定，且 beta 追蹤項嚴格凸，因此可行集合上的原問題仍有唯一 primal 解。

## 第 4 課：實際呼叫 MATLAB

從 MATLAB Current Folder 切到專案根目錄，輸入 `run('matlab/run_reconstruction.m')`。先看第 1–3 節如何建矩陣，第 4 節才求解。呼叫 `[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],opts);`，空陣列表示沒有一般不等式 A*x<=b；有上下限 lb/ub。x(1:121) 是 beta，x(122:241) 是 alpha。`sqrt(max(beta,0))` 僅用於顯示速度，原始 x 不可被任意裁剪後拿去做 KKT。

## 第 5 課：KKT 不是只看 exitflag

令等式乘子為 nu，下界與上界乘子為 lambdaL、lambdaU。Lagrangian = 0.5*x'*H*x+f'*x+nu'*(Aeq*x-beq)+lambdaL'*(lb-x)+lambdaU'*(x-ub)。對 x 微分得到 stationarity: H*x+f+Aeq'*nu-lambdaL+lambdaU=0。

四項 KKT：
1. Primal feasibility：Aeq*x=beq，且 lb<=x<=ub。
2. Dual feasibility：lambdaL,lambdaU>=0（nu 不受符號限制）。
3. Stationarity：上面的梯度為零。
4. Complementary slackness：lambdaL.*(x-lb)=0，lambdaU.*(ub-x)=0。

MATLAB 的 lambda.eqlin、lambda.lower、lambda.upper 對應這三組乘子。每項以最大絕對殘差或最大違反量檢查。這是凸 QP，KKT 成立是全域最優的充分條件。數值殘差只代表在容許誤差內成立。Active constraint 不一定有正乘子。

## 第 6 課：看懂四組結果

用 `readtable('results/independent/summary.csv')` 讀取已實際求解的獨立驗證結果。MATLAB 實際執行後請改讀 `results/matlab/summary.csv`。

ws=0 時加速度上限 48 段、煞車下限 36 段，共 84；ws=100 時分別 24、5，共 29。判斷方式是 `abs(alpha-1.5)<1e-4` 與 `abs(alpha+2)<1e-4`。wTrack=0.08 固定，因此 raw tracking 和 weighted tracking 的相對增幅相同，但數字不同。四個停車點就貢獻 2500 raw tracking，這是你面試時可以自己推導的 sanity check。

## 第 7 課：練習自己重寫（不改 V1.0 封存結果）

先在另一個臨時 .m 檔測試：
- 手寫 3 個 beta、2 個 alpha 的 Aeq 動力列；ds=0.5，第一列 [-1 1 0 -1 0]，第二列 [0 -1 1 0 -1]。
- 改成 4 個 alpha 手寫 D，驗證 D*alpha 是否等於 diff(alpha)。
- 說明 H 為什麼乘 2，以及為什麼 f 的 alpha 部分是 0。
- 用最小案例先求解，再回到完整 241 維程式。自行從 x 算出 E_track 和 E_smooth，對照 CSV，不要直接背結果。
- 說出此專案與路徑規劃的差異：幾何路徑已固定；本研究只安排沿路徑的縱向速度。它不是 Nav2 控制器，也不是實際購物車硬體測試。

## 面試表述（僅在自己能解釋後使用）

Original course project: Formulated a convex quadratic program for fixed-path vehicle speed planning with mandatory stops. Reconstruction: Rebuilt the documented MATLAB formulation, independently reproduced four weight settings using SciPy and checked full-space KKT residuals. Identified a weighted-versus-unweighted tracking metric labeling discrepancy. MATLAB execution verification is pending until personally run.

面試可以講研究取捨、模型限制、數據驗證與復現流程。不要宣稱最低行車時間、真正 jerk 最小化、實車驗證或已完成新的控制器。
