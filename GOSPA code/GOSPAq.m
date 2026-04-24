function [d_gospa, x_to_y_assignment, decomposed_cost] = ...
    GOSPAq(x_mat, y_mat, p, c, rho)
% This function implements the GOSPA quasimetric (q-metric) and is a modification of
% the GOSPA metric implementation by Abu Sajana Rahmathullah.
% The GOSPA q-metric generalises GOSPA alpha=2, and has a parameter rho

%A. F. García-Fernández, J. Gu, L. Svensson, Y. Xia, J. Krejčí, O. Kost, O. Straka 
% “GOSPA and T-GOSPA quasi-metrics for evaluation of multi-object tracking algorithms,” 
% accepted in IEEE Transactions on Aerospace and Electronic Systems, 2026.

%
%  [d_gospa, x_to_y_assignment] = GOSPA(x_mat, y_mat, p, c, alpha)
% computes the generalized optimal sub-pattern assignment metric (GOSPA)
% metric between the two finite sets, x_mat and y_mat for the given
% parameters c, p and alpha. Note that this implementation is based on
% auction algortihm, implementation of which is also available in this
% repository. For details about the metric, check
% https://arxiv.org/abs/1601.05585. Also visit https://youtu.be/M79GTTytvCM
% for a 15-min presentation about the metric.
%
% INPUT:
%   x_mat, y_mat: Input sets represented as real matrices, where the
%                 columns represent the vectors in the set
%   p           : 1<=p<infty, exponent
%   c           : c>0, cutoff distance
%   rho       : 0<rho<=1:  fraction of the maximum localisation error c^{p} that represents a false object cost (to the p-th power).
%
% OUTPUT:
%   d_gospa          : Scalar, GOSPA distance between x_mat and y_mat
%   x_to_y_assignment: Integer vector, of length same as the number of
%                      columns of x_mat. The i^th entry in the vector
%                      denotes the column index of the vector in y_mat,
%                      that is assigned to the vector in the i^th column of
%                      x_mat. Note that these indices are based on the
%                      permutation. Therefore, if #columns in x_mat <=
%                      #columns in y_mat, the entries in this vector will
%                      be between 1 and #columns in y_mat. Otherwise, the
%                      entries will be between 0 and #columns in y_mat,
%                      where 0 indicates that the corresponding columns in
%                      x_mat are unassigned.
%   decomposed_cost : Struct that returns the decomposition of the GOSPA
%                     metric for alpha=2 into 3 components:
%                          'localisation', 'missed', 'false'.
%                     Note that
%                     d_gospa = (decomposed_cost.localisation +
%                                decomposed_cost.missed       +
%                                decomposed_cost.false)^(1/p)
%
% Note: Euclidean base distance between the vectors in x_mat and y_mat is
% used in this function. One can change the function 'computeBaseDistance'
% in this function for other choices.

% check that the input parameters are within the valid range


n_ouput_arg=nargout;

checkInput();

nx = size(x_mat, 2); % no of points in x_mat
ny = size(y_mat, 2); % no of points in y_mat

% compute cost matrix
cost_mat = zeros(nx, ny);
for ix = 1:nx
    for iy = 1:ny
        cost_mat(ix, iy) ...
            = min(computeBaseDistance(x_mat(:, ix), y_mat(:, iy)), c);
    end
end

% intialise output values
decomposed_cost     = struct( ...
    'localisation', 0, ...
    'missed',       0, ...
    'false',        0);


x_to_y_assignment   = [];
opt_cost            = 0;

%Cost for false and a missed target
cost_false = rho*c^p;
cost_missed = (1-rho)*c^p;



% below, cost is negated to make it compatible with the 2-D asssignment algorithm
if nx == 0 % when x_mat is empty, all entries in y_mat are false
    opt_cost              = -ny * cost_false;
    decomposed_cost.false = opt_cost;
else
    if ny == 0 % when y_mat is empty, all entries in x_mat are missed
        opt_cost               = -nx * cost_missed;
        decomposed_cost.missed = opt_cost;

    else % when both x_mat and y_mat are non-empty, use 2-D assignment algorithm
        cost_mat = -(cost_mat.^p);

        [x_to_y_assignment, y_to_x_assignment]  = assign2D(cost_mat, true);


        % use the assignments to compute the cost
        for ind = 1:nx
            if x_to_y_assignment(ind) ~= 0
                opt_cost = opt_cost + cost_mat(ind,x_to_y_assignment(ind));



                decomposed_cost.localisation = ...
                    decomposed_cost.localisation ...
                    + cost_mat(ind,x_to_y_assignment(ind)) ...
                    .* double(cost_mat(ind,x_to_y_assignment(ind)) > -c^p);

                decomposed_cost.missed      = decomposed_cost.missed ...
                    - cost_missed ...
                    .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);

                decomposed_cost.false       = ...
                    decomposed_cost.false ...
                    - cost_false ...
                    .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);

            else
                opt_cost = opt_cost - cost_missed;

                decomposed_cost.missed = decomposed_cost.missed - cost_missed;

            end
        end
        opt_cost = opt_cost - sum(y_to_x_assignment == 0) * cost_false;

        decomposed_cost.false = decomposed_cost.false ...
            - sum(y_to_x_assignment == 0) * cost_false;

    end
end

% final output
d_gospa                      = (-opt_cost)^(1/p);
decomposed_cost.localisation = (-decomposed_cost.localisation);
decomposed_cost.missed       = (-decomposed_cost.missed);
decomposed_cost.false        = (-decomposed_cost.false);

    function checkInput()
        if size(x_mat, 1) ~= size(y_mat, 1)
            error('The number of rows in x_mat & y_mat should be equal.');
        end
        if ~((p >= 1) && (p < inf))
            error('The value of exponent p should be within [1,inf).');
        end
        if ~(c>0)
            error('The value of base distance c should be larger than 0.');
        end

        if ~((rho > 0) && (rho < 1))
            error('The value of rho should be within (0,1).');
        end
    end
end

function db = computeBaseDistance(x_vec, y_vec)
db = sum((x_vec - y_vec).^2)^(1/2);
end