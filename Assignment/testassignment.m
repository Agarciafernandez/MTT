function [samecost, sameassignment] = testassignment
%TESTASSIGNMENT  Test and compare assignment algorithms.
%   [SAMECOST, SAMEASSIGNMENT] = TESTASSIGN randomly generates distance
%   matrices and solves the assignment problem using different algorithms.
%   Edit the header of this file to change the simulation parameters.
%
%   <a href="assignment.html">assignment.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/6543">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=EVW2A4G2HBVAU">Donate via PayPal</a>
%
%   Markus Buehren
%   Last modified 14.03.2014
%
%   See also ASSIGNMENTOPTIMAL, ASSIGNMENTSUBOPTIMAL1,
%   ASSIGNMENTSUBOPTIMAL2, ASSIGNMENTALLPOSSIBLE.

% maximum simulation time in seconds
testTime = 20;       

% maximum matrix dimensions
minOrders = [ 1,  1];
maxOrders = [15, 25];

% If infAllowed in set to false, only distance matrices without infinite 
% costs are used
infAllowed = true;     

% If the product of the dimensions of the distance matrix is smaller than
% maxDimProduct, the assignment function computing all possible assignments
% (brute force solver) is used as reference. Set this to inf to use always or to
% 0 to never use.
maxDimProduct = 50; 

% use profiler or not
useProfiler = true;


% set recurstion limit
recursionLimit = get(0,'RecursionLimit');
set(0, 'RecursionLimit', 1000);

% start profiler
if useProfiler
  profile clear
  profile on
end

% initialize
startTime = clock;
nassign = 0;
h = waitbar(0, 'Please wait', 'Name', mfilename);

while 1
  for dim1 = minOrders(1):maxOrders(1)
    for dim2 = minOrders(2):maxOrders(2)
      if ~infAllowed
        
        % generate distMatrix without infinite elements
        distMatrix = rand(dim1,dim2);
        
      else
        
        if rand(1) < 0.5
          % generate distMatrix with some infinite elements
          distMatrix = rand(dim1,dim2);
          
          if rand(1) < 0.5
            % set some elements to inf
            distMatrix(rand(dim1,dim2) > rand(1)) = inf;
          end
          
        else
          
          % generate distMatrix with many infinite elements
          distMatrix = repmat(inf, dim1, dim2);
          
          if rand(1) < 0.5
            for row = 1:dim1
              if rand(1) < 0.8
                % set one element per row to finite number
                distMatrix(row, 1+floor(dim2*rand(1))) = rand(1);
              end
              if rand(1) < 0.3
                % set another element per row to finite number
                distMatrix(row, 1+floor(dim2*rand(1))) = rand(1);
              end
            end
          else
            for col = 1:dim2
              if rand(1) < 0.8
                % set one element per column to finite number
                distMatrix(1 + floor(dim1 * rand(1)), col) = rand(1);
              end
              if rand(1) < 0.3
                % set another element per column to finite number
                distMatrix(1 + floor(dim1 * rand(1)), col) = rand(1);
              end
            end
          end
          
        end
      end
      
      % quantize distMatrix
      if rand(1) < 0.3
        distMatrix = round(max(10, 1000*rand(1)^2)*distMatrix);
      end
      
      % transpose distMatrix
      if rand(1) < 0.5
        distMatrix = distMatrix';
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [assignmentCell{1}, costCell{1}] = assignmentoptimal    (distMatrix); %#ok
      [assignmentCell{2}, costCell{2}] = assignmentsuboptimal1(distMatrix); %#ok
      [assignmentCell{3}, costCell{3}] = assignmentsuboptimal2(distMatrix); %#ok
      if exist('munkres.m', 'file')
        [assignmentCell{4}, costCell{4}] = munkres_wrap(distMatrix); %#ok
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if dim1 * dim2 < maxDimProduct
        % if dimensions are moderate, compute all possible solutions
        [assignment_ref, cost_ref] = assignmentallpossible(distMatrix);
      else
        % use first algorithm as reference
        assignment_ref = assignmentCell{1};
        cost_ref       = costCell{1};
      end
      
      % count trial runs
      nassign = nassign + 1;
      
      if ~exist('samecost', 'var')
        M = length(costCell);
        samecost       = zeros(M,1);
        sameassignment = zeros(M,1);
      end
      
      % compute penalty for non-assignments
      finiteIndex = isfinite(distMatrix);
      penalty = max(max(distMatrix(finiteIndex))) * dim1 * dim2;
      cost_ref_mod = cost_ref;
      if ~isempty(penalty)
        cost_ref_mod = cost_ref + length(find(~assignment_ref)) * penalty;
      end
      
      costCell_mod = costCell;
      for m = 1:M
        
        % check validity of computed assignment
        if all(size(assignmentCell{m}) == [size(distMatrix, 1), 1]) && ...
            all(assignmentCell{m} <= size(distMatrix, 2)) && ...
            all(assignmentCell{m} >= 0)
          
          % check computed cost
          cost = 0;
          for k = 1:size(distMatrix, 1)
            if assignmentCell{m}(k) > 0
              cost = cost + distMatrix(k, assignmentCell{m}(k));
            end
          end
          
          if abs(costCell{m} - cost) < 100*eps
            
            % penalize non-assignments
            if ~isempty(penalty)
              costCell_mod{m} = costCell_mod{m} + length(find(~assignmentCell{m})) * penalty; %#ok
            end
            
            % compare costs
            if costCell_mod{m} <= cost_ref_mod + 100*eps
              samecost(m) = samecost(m) + 1;
            end
            
            % compare assignment with reference
            if all(assignmentCell{m} == assignment_ref)
              sameassignment(m) = sameassignment(m) + 1;
            end

          end % if abs(costCell{m} - cost) < 100*eps
        end % if all(size(assignmentCell{m}) == [dim1, 1])
      end % for m=1:M
    end % for dim2 = minOrder(2):maxOrders(2)
    
    % set waitbar
    curTime = etime(clock, startTime);
    waitbar(curTime/testTime, h);
  end
  
  % stop after given time
  if curTime > testTime
    break
  end
end
close(h);

% scale counters
samecost       = samecost/nassign;
sameassignment = sameassignment/nassign;

if useProfiler
  profile off
  profile report
end

% reset recurstion limit
set(0, 'RecursionLimit', recursionLimit);
