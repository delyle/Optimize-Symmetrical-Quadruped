function B = inpaint_nans_const(A)
% This function fills in NaNs that are surrounded by constant values.
% Useful for 2D matrices of categorical data stored as integers
% 
%%% Algorithm:
% 1: fill in any values that are surrounded by the same number or NaN
% 2: once (1) completed, fill in any NaN values that are surrounded by one
%    kind of number with one exception allowed
%
%%% Example:
% A = ...
%     [4 NaN NaN 4 4 4 2 NaN 2;
%      4 NaN NaN 4 4 NaN 2 2 2;
%      4 4 4 4 3 2 2 1 1;
%      3 3 3 3 4 1 1 NaN 1;
%      3 NaN 3 3 1 NaN 1 NaN 1;
%      3 NaN 3 4 1 NaN 1 NaN 1];
% B = inpaint_nans_const(A)
%
%%% output:
%
% B =
% 
%      4     4     4     4     4     4     2     2     2
%      4     4     4     4     4   NaN     2     2     2
%      4     4     4     4     3     2     2     1     1
%      3     3     3     3     4     1     1     1     1
%      3     3     3     3     1     1     1     1     1
%      3     3     3     4     1     1     1     1     1
% 
% Note that all values are filled, except for one, where it is not clear
% whether the value is a 4, 2 or 3.
%


% 1: fill in any values that are surrounded by the same number or NaN

[R,C] = find(isnan(A)); % Find NaN values
cNaN = NaN(size(A,1),1); % make column of NaNs
A_exp = [cNaN,A,cNaN]; % border A with columns of nans
rNaN = NaN(1,size(A_exp,2)); % make row of NaNs
A_exp = [rNaN;A_exp;rNaN]; % Border A (top, bottom) with rows of nans

for i = 1:length(R)
   % we add a row and column index to account for new NaN border
   r = R(i) + 1; c = C(i) + 1; 
   
   v = A_exp(r-1:r+1,c-1:c+1); % get surrounding values
   v = v(~isnan(v)); % ignore nans
   
   if ~isempty(v) && all(v == v(1))
       % all the values are the same
       A_exp(r,c) = v(1); % fill in the missing data with surrounding value
   end
end

% 2: fill in any values that are surrounded by the same number with one 
%    exception allowed

[R,C] = find(isnan(A_exp(2:end-1,2:end-1)));

for i = 1:length(R)
   % we add a row and column index to account for new NaN border
   r = R(i) + 1; c = C(i) + 1; 
   
   v = A_exp(r-1:r+1,c-1:c+1); % get surrounding values
   v = v([1:4,6:9]); % ignore middle value
   [vu,~,i_vu] = unique(v);
   vu_counts = accumarray(i_vu,1);
   
   if length(vu_counts) < 3 && min(vu_counts) == 1
       % there are no more than two unique values, and the smallest unique
       % value has only one occurance.
       % Now, fill in the missing data with the surrounding value that
       % occurs more than once
       A_exp(r,c) = vu(vu_counts > 1); 
   end
end

% Return B
B = A_exp(2:end-1,2:end-1);