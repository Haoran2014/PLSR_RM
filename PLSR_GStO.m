function [W,U]= PLSR_GStO(X0,Y0,p)
n = size(X0,2);
beta =5e-0; % beta can be adjusted.
Cx = X0'*X0+beta*eye(n);

Cxy = X0'*Y0;
Cy = Y0'*Y0;
rand('state',0)
W = orth(rand(n,p))/50;
U = orth((rand(size(Y0,2),p)))/2;
problem1.M = stiefelgeneralizedfactory(n, p, Cx);
problem1.cost = @cost;
    function f = cost(W)
        f = -trace(W'*Cxy*U);
    end
problem1.egrad = @egrad;
    function G = egrad(W)
        G = -Cxy*U;
    end

problem2.M = obliquefactory(size(Y0,2),p);
problem2.cost = @cost2;
    function f = cost2(U)
        f = -trace(W'*Cxy*U);
    end
problem2.egrad = @egrad2;
    function G = egrad2(U)
        G = -Cxy'*W;
    end

for i =1:10
    %[Xcg, xcost, info, options] = conjugategradient(problemM, X0);
    % checkgradient(problem1);    checkgradient(problem2);pause;
    
    
    [W, costw, info1, options1] = conjugategradient(problem1, W);
    [U, costu, info2, options2] = conjugategradient(problem2, U);
    
end
end