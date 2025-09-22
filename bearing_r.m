
p = [0 0; 1 0; 1 1; 0 1]; edges = [1 2; 2 3; 3 4; 4 1; 1 3]; % expected: bearing-rigid
p = [0 0; 1 0; 2 0; 3 0]; edges = [1 2; 2 3; 3 4]; % expected: not bearing-rigid
p = [0 0; 1 0; 0 1; -1 0; 0 -1]; edges = [1 2; 1 3; 1 4; 1 5]; % expected: not bearing-rigid
p = [0 0; 1 0; 0.5 0.9; -0.4 0.8; -0.8 -0.3; 0.6 -0.5]; edges = [1 2;1 3;1 4;1 5;1 6;2 3;2 4;2 5;2 6;3 4;3 5;3 6;4 5;4 6;5 6];
p = [0 0; 1 0; 1 1; 0 1]; edges = [1 2; 2 3; 3 4; 4 1; 2 4]; % expected: not bearing-rigid
p = [0 0; 2 0; 2 2; 0 2;0.5 1; 1.5 1]; edges = [1 2; 2 3; 3 4; 4 1; 1 5; 4 5; 5 6; 2 6; 3 6]; % expected: not bearing-rigid
d = 2;

isRigid = checkBearingRigidity(p, edges, d);


function isRigid = checkBearingRigidity(p, edges, d)
% 判断网络 (G,p) 是否是 Infinitesimally Bearing Rigid
%
% 输入：
%   p     : 节点位置向量 (n x d)
%   edges : 边集合 (m x 2)，每行 [i,j] 表示一条边
%   d     : 空间维度 (2 或 3)
%
% 输出：
%   isRigid : 布尔值，true 表示网络是 Infinitesimally Bearing Rigid

    n = size(p,1);   % 节点数
    m = size(edges,1); % 边数
    
    % 构造 incidence matrix H (m x n)
    H = zeros(m, n);
    for k = 1:m
        i = edges(k,1);
        j = edges(k,2);
        H(k,i) =  1;
        H(k,j) = -1;
    end
    H
    % 计算 edge vectors 和 gk
    E = (H * p) % m x d，每行是 e_k
    g = E ./ vecnorm(E,2,2); % 方向向量 gk
    
    % 构造对角块矩阵 diag(Pg/||ek||)
    blocks = cell(m,1);
    for k = 1:m
        gk = g(k,:)';
        ek = E(k,:)';
        P = eye(d) - gk*gk';
        blocks{k} = (P / norm(ek));
    end
    D = blkdiag(blocks{:});
    
    % Bearing rigidity matrix RB
    RB = D * kron(H, eye(d));
    
    % 求 null space 维度
    nullDim = size(null(RB),2);
    
    % 理论上，trivial motions 维度 = d(平移) + 1(缩放)
    trivialDim = d + 1;
    
    % 判断是否 rigid
    isRigid = (nullDim == trivialDim)
end







