function [assignmentVector, cost] = munkres_wrap(distMatrix)

[assignment, cost] = munkres(distMatrix); %#ok
M = size(assignment, 1);
if all(size(assignment) == size(distMatrix))
  % wrapper for munkres.m from http://www.mathworks.com/matlabcentral/
  % fileexchange/20328-munkres-assignment-algorithm
  assignmentVector = zeros(M, 1);
  for k = 1:M
    index = find(assignment(k,:));
    if length(index) > 1
      keyboard
    end
    assert(length(index) <= 1);
    if ~isempty(index)
      assignmentVector(k) = index(1);
    end
  end
elseif all(size(assignment) == [1, size(distMatrix,1)])
  % wrapper from munkres.m from http://www.mathworks.com/matlabcentral/
  % fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3
  assignmentVector = assignment';
else
  % unknown version, return empty assignment
  assignmentVector = zeros(size(distMatrix,1), 1);
end