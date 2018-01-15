function [ T, X, sizeBlocks ] = jordanForm( A, type )
%jordanForm Gives the matrix transformation to the jordan from of A with
%the size of each block. 
%
%   OUTPUT
%   - T is the transformation matrix such that inv(T)*A*T is the real jordan
%   form of A.
%   - X is a sdp variable (using Yalmip). This is a matrix the size of A with
%   k diagonal blocks. k is the number of different eigenvalues of A. The
%   ith diagonal block of X has the same size that the algebraic
%   multiplicity of the corresponding eigenvalue in A.
%   - sizeBlocks is a list of the albegraic multiplicity of each eigenvalue
%   of A.
%
%   For example, if A = [0.2 0; 0.2 -0.2], then T = [2 0; 1 1] such that
%   inv(T)*A*T = [0.2 0; 0 -0.2]. Then sizeBlocks = [1, 1], because the
%   algebraic multiplicity of 0.2 is 1 and so does for -0.2. Then X is a
%   diagonal matrix. For more example, please refer to the article cited
%   below.
%
%   [T, X, sizeBlocks] = jordanForm(A) computes T, X and sizeBlocs as
%   defined previously.
%   [T, X, sizeBlocks ] = jordanForm(A, 'sof') computes T, X and sizeBlocks
%   as previously but X is diagonal and not block diagonal. It can be used
%   for static output feedback design for example. See the following
%   article for more information.
%
%   Version 1.0 / January 2018
% 
%   If you are using or modifying this code, please cite the following
%   reference:
%   M. Barreau, F. Gouaisbaut and A. Seuret,
%   Static tatic State and Output Feedback Synthesis for Time-Delay Systems
%
%   See also generateEpsilon

if nargin == 1
    sizeSdpVar = size(A, 1);
else
    if strcmpi(type, 'sof')
        sizeSdpVar = 1;
    else
        sizeSdpVar = size(A, 1);
    end
end

%% Jordan form of A
[T, A] = jordan(A); % Complex jordan form of A
[T, A] = cdf2rdf(T, A); % Real jordan form of A and T is the transition matrix

%% Creation of sizeBlocks
sizeBlocks = [];
lambda = A(1,1); % First eigenvalue of A
j = 0;
for i=1:size(A,1)
    if A(i,i) ~= lambda % We compare, if it is the same eigenvalue than before, the block continue
        sizeBlocks = [sizeBlocks, j]; % Size of the block
        lambda = A(i, i); % New eigenvalue
        j = 0;
    end
    j = j+1;
end
sizeBlocks = [sizeBlocks, j];

%% Creation of the X matrix
X = [];
for i=1:length(sizeBlocks)
    if sizeSdpVar > 1
        X = blkdiag(X, sdpvar(sizeBlocks(i), sizeBlocks(i), 'full'));
    else
        X = blkdiag(X, sdpvar(1)*eye(sizeBlocks(i)));
    end
end

end
