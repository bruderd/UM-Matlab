function w = get_w_from_beta( beta , horizon )
%get_w_from_beta: Converts Volterra series coefficients from regression
%                  to the format from the Isidori paper.
%   
% Inputs:
%   Alpha - Volterra coefficients from regression
%   horizon - discrete time horizon over which inputs are defined 

% dreallocate for sequence of w_j's, which will be stored as rows of w
w = zeros( horizon , 2^(horizon-1) );

% determine how many of each degree will be in each w_j
pascals_triangle = fliplr( pascal( horizon ) );
num_degree = zeros( horizon , horizon );
for i = 1 : horizon
    num_degree(i,1:i) = diag( pascals_triangle , -(i-1) + (horizon-1) )';
end

% allocate elements of Alpha in w
Bindex = 1; % keeps track of which element we are currently assigning
wassigned = zeros( size(w) ); % keeps track of which elements of w have been assigned

for i = 1 : horizon   % 1 sorting round per degree
    
    for bucket = horizon : -1 : 1   % start from the largest bucket
        num_assigned = 0;   % keep track of how many elements of Alpha have been assinged to the current bucket
        
        % choose the next empty spot in the bucket that follows the 1-2-1 rule
        while num_assigned < num_degree( bucket , i )
            for j = 1 : size( w , 2 ) % follow the 1-2-1 pattern
                if wassigned( bucket , j ) == 0 
                    % Case 0: It's in the very first spot
                    if j == 1
                        windex = j;
                        break;  % should stop if from goint into case 1
                    end
                    % Case 1: It's in a first spot
                    if mod(j,4) == 1 && wassigned(bucket,j-1) ~= i
                        windex = j;
                        break;
                    end
                    % Case 2: It's in a second spot
                    if mod(j,4) == 2 && wassigned(bucket,j-1) ~= i
                        windex = j;
                        break;
                    end
                    % Case 2: It's in a third spot
                    if mod(j,4) == 3
                        windex = j;
                        break;
                    end
                    % Case 4: It's in a fourth spot
                    if mod(j,4) == 0 && wassigned(bucket,j-1) ~= i
                        windex = j;
                        break;
                    end
                end
            end
            
            w( bucket , windex ) = beta( Bindex ); % assign value
            
            Bindex = Bindex + 1; % increment assigment counter
            wassigned( bucket , windex ) = i;   % mark that spot as assigned with coefficent for degree i term
            num_assigned = num_assigned + 1;
        end
    end
    
end


end

