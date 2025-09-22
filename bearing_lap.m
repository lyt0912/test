

edges = [1 2;1 3;1 4;1 5;1 6;2 3;2 4;2 5;2 6;3 4;3 5;3 6;4 5;4 6;5 6];
E = [-1.0000 0; -0.5000 -0.9000; 0.4000 -0.8000; 0.8000 0.3000; -0.6000 0.5000; 0.5000 -0.9000; 1.4000 -0.8000; 1.8000 0.3000; 0.4000 0.5000; 0.9000 0.1000; 1.3000 1.2000; -0.1000 1.4000; 0.4000 1.1000; -1.0000 1.3000; -1.4000 0.2000];
% 边连接
p = 6; %节点个数
d = 2;

isRigid = checkBearingRigidity(p, E, edges, d);


function isRigid = checkBearingRigidity(p, E, edges, d)

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
    H_ = kron(H, eye(d));
    
    g = E ./ vecnorm(E,2,2); % 方向向量 gk
    % 构造对角块矩阵 diag(Pg/||ek||)
    blocks = cell(m,1);
    for k = 1:m
        gk = g(k,:)';
        P = eye(d) - gk*gk';
        blocks{k} = (P );
    end
    D = blkdiag(blocks{:})
    
    B =  H_' * D' * D * H_;
    
    % 求 null space 维度
    nullDim = size(null(B),2);
    
    % 理论上，trivial motions 维度 = d(平移) + 1(缩放)
    trivialDim = d + 1;
    
    % 判断是否 rigid
    isRigid = (nullDim == trivialDim)
end