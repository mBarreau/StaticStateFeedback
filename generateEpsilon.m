function [ epsilon, subA ] = generateEpsilon( A, sizeBlocks )
%generateEpsilon Generates the E matrix corresponding to matrix A.
%   This function generates a matrix E such that the following hold:
%   - E is block diagonal and has as many diagonal block as
%   length(sizeBlocs);
%   - Each diagonal block of E has as many rows has its corresponding value
%   in sizeBlocks;
%   - Each diagonal block of E is lambda * I (with I the identity matrix of
%   the correct size) where lambda is the trace of the same block in A.
%
%   OUTPUT:
%   epsilon is the E matrix;
%   subA is a cell containing length(sizeBlocks) matrices which are the
%   corresponding diagonal blocks of A. See the code for more details.
%
%   [epsilon, subA] = generateEpsilon( A, sizeBlocks).
%   [epsilon, subA] = generateEpsilon( A ) produces the same than
%   generateEpsilon( A, length(A)).
%
%   Version 1.0 / January 2018
% 
%   If you are using or modifying this code, please cite the following
%   reference:
%   M. Barreau, F. Gouaisbaut and A. Seuret,
%   Static tatic State and Output Feedback Synthesis for Time-Delay Systems
%
%   See also jordanForm

switch nargin
    case 1
        sizeBlocks = size(A,1);
    case 2
        if sizeBlocks <= 0
            sizeBlocks = size(A,1);
        end
    otherwise
end

if size(A,1) ~= size(A, 2)
   error(' in generateEpsilon: A is not square.'); 
end

if sum(sum(sizeBlocks)) ~= length(A)
   warning(' in generateEpsilon: the sum of all values in sizeBlocks is not equal to length(A)');
end

m = length(sizeBlocks); % Number of diagonal blocks
subA = cell(1, m); % Cell containing the subBlocks of A
epsilon = zeros(size(A,1)); % The E matrix

for i=1:m
    indexStart = sum(sizeBlocks(1:(i-1)))+1;
    indexEnd = indexStart + sizeBlocks(i) -1;
    range = indexStart:indexEnd; % This is the range of a diagonal block of A
    subA{i} = A(range, range); % Extraction of the block
    epsilon(range, range) = trace(subA{i})*eye(sizeBlocks(i))/sizeBlocks(i); % Creation of the E matrix
end

end

