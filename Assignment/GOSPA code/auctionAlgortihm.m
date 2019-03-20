function [person_to_obj, obj_to_person, opt_cost] ...
    = auctionAlgortihm(cost_mat, max_iter)
% AUTHOR: Abu Sajana Rahmathullah
% DATE OF CREATION: 7 August, 2017
%
%  [person_to_obj, obj_to_person, opt_cost] ...
%                   = auctionAlgortihm(cost_mat, epsil, max_iter)
% returns the results of auction algorithm in assigning 'm' unique persons
% to the 'n' objects, such that the overall cost is maximised.
%
% REFERENCE: Section 6.5.1 in 'Design and Analysis of Modern Tracking
%             Systems' by Blackman and Popoli, 1999 edition
%
% INPUT:
%   cost_mat: matrix of size mxn consising of cost of assigning one of m
%             persons to one of each n objects
%   max_iter: maximum number of iterations to terminate the alg
%
% OUTPUT:
%   person_to_obj: A vector of size 1xm which has indices of the assigned
%                  objects or '0' if unassigned
%   obj_to_person: A vector of size 1xn, which has indices of the
%                  assigned persons and a vector
% 	opt_cost     : A scalar, value of the (optimal) assignment made
%

m_person = size(cost_mat, 1); n_obj = size(cost_mat,2);

% index for objects that will be left unassigned
un_ass_ind = 0;

% corner cases when there is only one person and/or one object
if m_person == 1
    [opt_cost, person_to_obj] = max(cost_mat);
    obj_to_person = un_ass_ind * ones(n_obj, 1);
    obj_to_person(person_to_obj) = 1;
    return;
end
if n_obj == 1
    [opt_cost, obj_to_person] = max(cost_mat);
    person_to_obj = un_ass_ind * ones(m_person, 1);
    person_to_obj(obj_to_person) = 1;
    return;
end

swap_dim_flag = false;
epsil = 1/ max(n_obj, m_person);

if n_obj < m_person % the below implementation works for m_person<=n_obj.
    % If not satisifed, swap the cost matrix
    cost_mat = cost_mat';
    m_person = size(cost_mat, 1); n_obj = size(cost_mat,2);
    swap_dim_flag = true;
end


% obj ind of assignment for person 1 to m
person_to_obj = un_ass_ind * ones(m_person, 1);

% person ind of assignemnt for obj 1 to n
obj_to_person = un_ass_ind * ones(n_obj, 1);
opt_cost = 0;
p_obj = zeros(1, n_obj); % price for each object

iter = 0;
while(~all(person_to_obj ~= un_ass_ind))
    if iter > max_iter
        warning(['Maximum number of iterations reached! Retry with ' ...
            'more than ' num2str(max_iter) ' number of iterations.']);
        break;
    end
    for i = 1:m_person
        if(person_to_obj(i) == un_ass_ind)
            % pick the unassigned person i to bid for its best obj jStar
            
            % value for each object for the person i
            [val_i_j, j] = sort(cost_mat(i, :) - p_obj, 2, 'descend');
            
            % j_star is the best object for person i
            j_star = j(1);
            
            % 1st and 2nd best value for person i
            v_ijStar = val_i_j(1);
            w_ijStar = val_i_j(2);
            
            
            % bid for obj jStar
            if w_ijStar ~= -Inf % if there is no 2nd best
                p_obj(j_star) = p_obj(j_star) + v_ijStar - w_ijStar + epsil;
            else
                p_obj(j_star) = p_obj(j_star) + v_ijStar + epsil;
            end
            
            
            if obj_to_person(j_star) ~= un_ass_ind % if j_star is uassigned
                % if j_star is assigned to obj_to_person(j_star) before,
                % remove the cost of old assignment
                % and unassign the person obj_to_person(j_star)
                opt_cost = opt_cost ...
                    - cost_mat(obj_to_person(j_star), j_star);
                person_to_obj(obj_to_person(j_star)) = un_ass_ind;
            end
            obj_to_person(j_star) = i; % assign i to j_star
            person_to_obj(i) = j_star; % assign j_star to i
            
            % update the cost of new assignemnt
            opt_cost = opt_cost + cost_mat(i, j_star);
        end
    end
    iter = iter+1;
end

% swap the dimensions to compensate for the swapping of the cost matrix
% done in the beginning
if swap_dim_flag == true
    tmp = obj_to_person;
    obj_to_person = person_to_obj;
    person_to_obj = tmp;
end
