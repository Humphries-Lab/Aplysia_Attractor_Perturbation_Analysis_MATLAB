function [R, D, epsOut] = fn_recurrencePlot(traj, epsilon, distType)
% FN_RECURRENCEPLOT  Minimal pairwise recurrence-plot matrix, for VISUAL
% display only (imagesc(R)). Does not affect any recurrence-density number
% computed elsewhere -- those all come from fn_recurrenceDensity /
% fn_crossRecurrenceDensity.
%
% INPUTS
%   traj     - N x D matrix, trajectory points
%   epsilon  - recurrence threshold
%   distType - 'euc' (only Euclidean is implemented here)
%
% OUTPUTS
%   R      - N x N logical recurrence matrix
%   D      - N x N distance matrix
%   epsOut - epsilon used (passthrough, for convenience)
 
if nargin < 3, distType = 'euc'; end
if ~strcmpi(distType, 'euc')
    warning('fn_recurrencePlot: only euclidean distance implemented, ignoring distType=%s', distType);
end
 
N = size(traj, 1);
D = zeros(N, N);
for i = 1:N
    D(i,:) = sqrt(sum((traj - traj(i,:)).^2, 2))';
end
R = D <= epsilon;
epsOut = epsilon;
 
end
