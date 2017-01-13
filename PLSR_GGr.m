function W= PLSR_GGr(X0,Y0,p)
n = size(X0,2);

beta =5e-0; % beta can be adjusted.
Cxx = X0'*X0+beta*eye(n); %5e-0

Cxy = X0'*Y0;
C = Cxy *Cxy';
%%
rand('state',0);

W = orth(rand(n,p))/100;
%%

    problem.M  = grassmanngeneralizedfactory(n, p,Cxx);
    problem.cost = @cost;
    function f = cost(W)
    f = -trace(W'*C*W);
    end
    problem.egrad = @egrad;
    function G = egrad(W)
    G = -2*C*W;
    end
%checkgradient(problem);% pause;
%checkhessian(problem); pause;
 

  [W, costw, info1, options1] = conjugategradient(problem, W);
    
    
  %[W, costw, info1, options1] = steepestdescent(problem, W);
    
  

end