function [d_gospa, x_to_y_assignment, decomposed_cost] = ...
    GOSPA(x_mat, y_mat, p, c, alpha)
% AUTHOR: Abu Sajana Rahmathullah
% DATE OF CREATION: 7 August, 2017
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
%   alpha       : 0<alpha<=2, factor for the cardinality penalty.
%                 Recommended value 2 => Penalty on missed & false targets
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

dummy_cost = (c^p) / alpha; % penalty for the cardinality mismatch

% below, cost is negated to make it compatible with auction algorithm
if nx == 0 % when x_mat is empty, all entries in y_mat are false
    opt_cost              = -ny * dummy_cost;
    decomposed_cost.false = opt_cost;
else
    if ny == 0 % when y_mat is empty, all entries in x_mat are missed
        opt_cost               = -nx * dummy_cost;
        
        if(alpha==2)
            decomposed_cost.missed = opt_cost;
        end
    else % when both x_mat and y_mat are non-empty, use auction algorithm
        cost_mat = -(cost_mat.^p);
        [x_to_y_assignment, y_to_x_assignment, ~] ...
            = auctionAlgortihm(cost_mat, 10*(nx * ny));
        % use the assignments to compute the cost
        for ind = 1:nx
            if x_to_y_assignment(ind) ~= 0
                opt_cost = opt_cost + cost_mat(ind,x_to_y_assignment(ind));
                
                if(alpha==2)
                    
                    decomposed_cost.localisation = ...
                        decomposed_cost.localisation ...
                        + cost_mat(ind,x_to_y_assignment(ind)) ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) > -c^p);
                    
                    decomposed_cost.missed      = decomposed_cost.missed ...
                        - dummy_cost ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);
                    
                    decomposed_cost.false       = ...
                        decomposed_cost.false ...
                        - dummy_cost ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);
                end
            else
                opt_cost               = opt_cost - dummy_cost;
                if(alpha==2)
                    decomposed_cost.missed = decomposed_cost.missed - dummy_cost;
                end
            end
        end
        opt_cost = opt_cost - sum(y_to_x_assignment == 0) * dummy_cost;
        if(alpha==2)
            decomposed_cost.false = decomposed_cost.false ...
                - sum(y_to_x_assignment == 0) * dummy_cost;
        end
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
        
        if ~((alpha > 0) && (alpha <= 2))
            error('The value of alpha should be within (0,2].');
        end
        if alpha ~= 2 && n_ouput_arg==3
            warning(['decomposed_cost is not valid for alpha = ' ...
                num2str(alpha)]);
        end
    end
end

function db = computeBaseDistance(x_vec, y_vec)
db = sum((x_vec - y_vec).^2)^(1/2);
end