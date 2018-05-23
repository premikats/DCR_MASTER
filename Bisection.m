function [x,iter,flag] = Bisection(fnhandle,x_l,x_u,varargin)
%Bisection
% Bisection method is used for finding the root of the given function

%%
% Inputs:
% x_l = lower bound
% x_u = upper bound

% Outputs:
% x = Root
% iter = no. of iterations to obtain the solution
% flag = status of the problem
%           0: solved
%           1: improper bounds
%           2: not solved within the given no. of iterations
% Last argument - no. of iteration
% Last but before - Tolerance

%%

% DON'T BE SAD PLEASEEEEEEEEEE
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n_f = nargin(fnhandle);
n = nargin;

if n<n_f+2||n>n_f+4
    error('number of arguments do not match');
end

tol=10*eps;
N=1000000;

if n==n_f+3
    tol = varargin{n_f};
elseif n==n_f+4
    tol = varargin{n_f};
    N = varargin{n_f+1};
end
count = 1;
flag = 2;
while N>0
    f_l = fnhandle(x_l,varargin{1:n_f-1});
    f_u = fnhandle(x_u,varargin{1:n_f-1});
    if f_l*f_u>0
        disp('The root does not lie within the lower and upper bounds given');
        x = 0;
        iter = count;
        flag = 1;
        return
    else
        c = (x_l+x_u)/2;
        if (x_u-x_l)<tol
            x = c;
            iter = count;
            flag = 0;
            return
        end
        f_c = fnhandle(c,varargin{1:n_f-1});
        if f_c*f_l<0
            x_u = c;
        elseif f_c*f_u<0
            x_l = c;
        elseif f_c==0
            x = c;
            iter = count;
            flag = 0;
            return
        end
    end
    count = count+1;
    N = N-1;
end
if flag==2
    disp('Unable to solve within max number of iterations');
    x = 0;
    iter = count;
end
        
        
