function [x,x_list,num_iters] = fossils_test(A,b,varargin)
%FOSSILS Solve A*x = b in the least-squares sense using FOSSILS
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 12*size(A,2))
%   - iterations: can be set in three ways
%        * if iterations is a vector of integers, iterations(i) contains
%          the number of steps to be taken on the i-th refinement step.
%        * if iterations is a non-integer, iterations sets the tolerance
%          for the backward error and the number of steps is determined
%          adaptively
%        * if iterations is 'adaptive', the number of steps is determined
%          adaptively with the default tolerance (default 'adaptive')

    tic
    
    Anum = isnumeric(A);
    if Anum
        scale = vecnorm(A);
        A = A ./ scale;
        m = size(A,1);
        n = size(A,2);
        Afun = @(x,op) mul(A,x,op);
    else
        scale = 1;
        m = size(b,1);
        n = size(A(zeros(size(b,1),0), true),1);
        Afun = A;
    end

    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = 12*n;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        iterations = varargin{2};
    else
        iterations = 'adaptive';
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3})
        S=varargin{3};
    end

   
    adaptive = true;
    iterations = [100,100];
    betol = 1e-20;


    if Anum
        SA = full(S*A);
    else
        SA = full(Afun(S',true)');
    end
    [U,svals,V] = svd(SA,'econ'); svals = diag(svals);
    x = V*((U'*(S*b)) ./ svals);
    x_list = x;
    r = b - Afun(x,false);

    % Chebyshev interpolation
    momentum = n/d;
    damping = (1 - momentum)^2;

    % Norm and condition number estimation
    Acond = max(svals) / min(svals);
    Anorm = max(svals); 
    Afronorm = norm(svals);
    bnorm = norm(b);

    if Acond > 1/eps/30
        reg = 10 * Afronorm * eps;
        sreg = sqrt(svals.^2 + reg^2);
    else
        reg = 0;
        sreg = svals;
    end

    betol = Afronorm * betol; % Backward error tolerance

    num_iters = 0;

    for loop = 1:length(iterations)
        c = (V'*(Afun(r,true) - reg^2*x))./sreg;
        dy = c; dyold = dy;
        for i = 1:iterations(loop)
            num_iters = num_iters + 1;
            z = V*(dy./sreg);
            update = damping * (c - (V'*(Afun(Afun(z,false),true)...
                + reg^2 * z)) ./sreg) + momentum*(dy-dyold);
            dyold = dy;
            dy = dy + update;

            if true
                xhat = x + V*(dy./sreg);
                rhat = b - Afun(xhat,false);
                be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
                fprintf('Iteration %d\t%e\t%e\n', i, norm(update), be);
            end
%             if adaptive && loop == 1 && norm(update) ...
%                     <= 10*(Anorm*norm(x) + 0.04*Acond*norm(r))*eps
%                 break
%             elseif adaptive && loop == 2 && mod(i,5) == 0
%                 xhat = x + V*(dy./sreg);
%                 rhat = b - Afun(xhat,false);
%                 be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
%                 if be <= betol; break; end
%             end
            xhat = x + V*(dy./sreg);
            x_list(:,end+1)=xhat;
        end
        x = x + V*(dy./sreg);
        r = b - Afun(x,false);
        if adaptive
            be = posterior_estimate(Afun,x,r,V,svals,Afronorm,bnorm);
            if be <= betol; break; end
            if loop == 2
                warning('Exiting without backward stability!')
            end
        end
    end
    x_list = x_list./ scale.';
    x = x ./ scale.';
end

function be = posterior_estimate(Afun,x,r,V,svals,Afronorm,bnorm)
    theta = Afronorm / bnorm;
    xnorm = norm(x);
    be = norm((V'*Afun(r,true)) ./ (svals.^2 ...
        + theta^2*norm(r)^2/(1+theta^2*xnorm^2)).^0.5) ...
        * (theta / sqrt(1+theta^2*xnorm^2));
end

function b = mul(A, x, op)
    if op
        b = A'*x;
    else
        b = A*x;
    end
end