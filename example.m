%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Exemple of solution for a BMI using       %
%          the jordan form methodology             %
%   15/01/2018              by M. Barreau          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This exemple is a simple implementation of the methodology developped in
% the conference paper: Static State and Output Feedback Synthesis for
% Time-Delay Systems. It solves a BMI and computes a static gain feedback
% using YALMIP with the sdpt3 solver. For more information about the LMI
% used, please refer to the cited article. The problem we solve has the
% following structure: x'*Phi0*x < 0 for all Mx = 0. This condition is
% ensured if the following LMI is satisfied: Phi0 + He(M'*FW) < 0. The aim
% of this short code is to find a "good" value of FW.

clear
close all
clc

%% System parameters
A = [0.2 0; 0.2 -0.2];
B = [-1 0; -1 -1];
K = [1.2 0; -1 1.8]; % Initial point
N = 1; % Order of the LMI
h = 1e-5; % Time-delay
lmax = 3; % Number of iterations
nbDecimal = 3; % Precision of the path-following algorithm

%% Creation of the BMI

% System transformation into its jordan form
[ T, X, sizeBlocs ] = jordanForm(A);
A = (T\A)*T;
B = T\B;
K = K*T;

% Solver parameters
n = length(A);
S = sdpvar(n);
R = sdpvar(n);
PN = sdpvar(n*(N+1));
Kbar = sdpvar(size(B, 2), n, 'full');

% Creation of matrices (only technical)
gammaN = @(N, k, i) -(2*i+1)*(1-(-1)^(k+i))*(i <= k);
GammaN = zeros((N+1)*n, n*(N+3));
for k = 0:N
    GammaN((k*n+1):((k+1)*n), 1:(3*n)) = [zeros(n), eye(n), (-1)^(k+1)*eye(n)];
    for i=0:(N-1)
        GammaN((k*n+1):((k+1)*n), ((3+i)*n+1):((4+i)*n)) = gammaN(N,k,i)*eye(n);
    end
end

Stilde = blkdiag(zeros(n), S, -S, zeros(n*N));
Rtilde = [];
for i=0:N
    Rtilde = blkdiag(Rtilde, (2*i+1)*R);
end

FN =  [eye(n) zeros(n, (N+2)*n)];
GN = [zeros(n) eye(n) zeros(n, (N+1)*n);
    zeros(n*N, 3*n), h*eye(n*N)];
HN = [FN', GammaN(1:n*N, :)']';
M = [eye(n)*X -A*X -B*Kbar zeros(n, n*N)];

He = GN'*PN*HN;
Phi0 = He + He' + Stilde + h^2*FN'*R*FN - GammaN'*Rtilde*GammaN;

%% Solving the BMI using an heuristic

disp('===================================');
disp(strcat(['Solving the BMI problem for h=', num2str(h),'.']));
disp('===================================');

option = sdpsettings('verbose', 0, 'solver', 'sdpt3');

for l = 1:lmax
    % Generation of the E matrices
    E1 = generateEpsilon(eye(n), sizeBlocs);
    E2 = generateEpsilon(-A, sizeBlocs);
    E3 = generateEpsilon(-B*K, sizeBlocs); % We generate using the initial guess
    FW = [E1 E2 E3 repmat(zeros(n)/n, 1, N)];
    
    T = M'*FW;
    T = T + T';
    
    Constraints = [PN >= 1e-5, S >= 1e-5, R >= 1e-5, Phi0+T <= -1e-5];
    
    optimize(Constraints, {}, option);
    pres = checkset(Constraints);
    
    if sum(pres > 0) < length(pres)
        disp(strcat(['No solution found at h=', num2str(h),'.']));
        break;
    else
        K = value(Kbar)/value(X); % Get the new gain K
    end
end

%% Example of path-following method
if sum(pres > 0) == length(pres)
    
    disp(strcat(['There is a solution found at h=', num2str(h),'.']));
    disp('===================================');
    disp('Starting the path-following search for hmax');
    disp('===================================');
    
    solved = -1;
    numPb = 0;
    epsilon = 1;
    h = h + epsilon;
    i = 0;
    hmax = 0;
    while i <= nbDecimal && h >= 0
        
        % New LMI conditions because h changed
        GN = [zeros(n) eye(n) zeros(n, (N+1)*n);
            zeros(n*N, 3*n), h*eye(n*N)];
        He = GN'*PN*HN;
        Phi0 = He + He' + Stilde + h^2*FN'*R*FN - GammaN'*Rtilde*GammaN;
        
        Constraints = [PN >= 1e-5, S >= 1e-5, R >= 1e-5, Phi0+T <= -1e-5];
        
        % Solve the problem
        optimize(Constraints, {}, option);
        pres = checkset(Constraints);
        
        disp(strcat(['----> Trying at h=', num2str(h),'.']));
        
        if(sum(pres > 0) < length(pres))
            solved = -1;
        else
            solved = 1;
            K = value(Kbar)/value(X);
            
            for l = 1:lmax
                % Update of the third E matrix
                E3 = generateEpsilon(-B*K, sizeBlocs); % We generate using the new value of K
                FW = [E1 E2 E3 repmat(zeros(n)/n, 1, N)]; % New FW matrix such that FW is close to M
                
                T = M'*FW;
                T = T + T';
                
                Constraints = [PN >= 1e-5, S >= 1e-5, R >= 1e-5, Phi0+T <= -1e-5];
                
                optimize(Constraints, {}, option);
                pres = checkset(Constraints);
                
                % There can be numerical problem and then sum(pres > 0) <
                % length(pres) even if theoritically it is not possible. Then
                % we count this as an error and if the error propagates, we
                % stop the algorithm
                if sum(pres > 0) < length(pres)
                    K = value(Kbar)/value(X);
                    numPb = numPb + 1;
                    if numPb > 1
                        solved = -1;
                        warning('Numerical problem. The algorithm has been avorted.');
                    end
                else % No numerical problem, we continue
                    numPb = 0;
                    K = value(Kbar)/value(X);
                end
            end
        end
        
        % This is the core of the path-following algorithm
	% We change the value of h depending on whether the problem was feasible or not
        if solved > 0
            disp('There is a solution. Increase h.');
            hmax = h;
            h = h + epsilon*10^(-i);
        else
            disp('No solution. Decrease h.');
            i = i+1;
            h = h+epsilon*(10^(-i)-10^(-i+1));
            if h <= 0
                h = 10^(-nbDecimal);
            end
        end
        
    end
    
    % End of the pth-following search
    disp('===================================');
    disp(strcat(['Maximum hmax=', num2str(hmax),'.']));
    disp('===================================');
end
