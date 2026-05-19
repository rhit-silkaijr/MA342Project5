function [Et,Mt,A] = Ising_2D_Func(T, A, sample, num_iters)
    % N = base lattice size, padded boundary included
    N = length(A);
    
    % steps
    t_max = num_iters;
    
    % random ass constants, idk
    J = 1;
    n = 1;
    
    % Neighbor mapping
    % square grid
    neighbors = [-1, 0; 1, 0; 0, 1; 0, -1]; % left, right, up, and down (dx, dy)

    if sample

        % num_iters_discarded
        t_pre = num_iters/5;
    
        for t = 1:t_pre
            A(1,:) = A(end-1,:);
            A(end,:) = A(2,:);
            A(:,1) = A(:,end-1);
            A(:,end) = A(:,2);
            % hard set corners
            A(1,1) = A(end-1,end-1);
            A(1,end) = A(end-1, 2);
            A(end,1) = A(2,end-1);
            A(end,end) = A(2,2);
        
            % get mask of % that can change
            change_mask = rand(size(A)) > 0.4;
            change_mask(1,:) = 0; change_mask(end,:) = 0;
            change_mask(:,1) = 0; change_mask(:,end) = 0;
        
            A_new = A;
            % loop through all indices that can change
            for x = 2:N-1
                for y = 2:N-1
                    if change_mask(y,x) == 1
                        % calculate E from surrounding neighbors (p. 422)
                        E = 0;
                        for i = 1:length(neighbors)
                            E = E + A(y+neighbors(i,2),x+neighbors(i,1))*A(y,x);
                        end
                        E = E*(-J);
                        % calculate E from neighbors as is A(y,x) was flipped
                        A_temp = A; A_temp(y,x) = -A_temp(y,x);
                        E_flipped = 0;
                        for i = 1:length(neighbors)
                            E_flipped = E_flipped + A_temp(y+neighbors(i,2),x+neighbors(i,1))*A_temp(y,x);
                        end
                        E_flipped = E_flipped*(-J);
                        % use equation 12.13 to decide if we keep the flip
                        if E > E_flipped
                            P = 1;
                        else
                            P = exp((E-E_flipped)/T);
                        end
                        if P > rand()
                            A_new(y,x) = -A(y,x);
                        end
                    end
                end
            end
            % set lattice for next iteration
            A = A_new;
        end
    end
    
    Et = zeros(t_max, 1);
    Mt = zeros(t_max, 1);
    for t = 1:t_max
        % hard set edges of A for torroidal boundary condition (square lattice)
        A(1,:) = A(end-1,:);
        A(end,:) = A(2,:);
        A(:,1) = A(:,end-1);
        A(:,end) = A(:,2);
        % hard set corners
        A(1,1) = A(end-1,end-1);
        A(1,end) = A(end-1, 2);
        A(end,1) = A(2,end-1);
        A(end,end) = A(2,2);
    
        % get mask of % that can change
        change_mask = rand(size(A)) > 0.4;
        change_mask(1,:) = 0; change_mask(end,:) = 0;
        change_mask(:,1) = 0; change_mask(:,end) = 0;
    
        A_new = A;
        E_tot = 0;
        % loop through all indices that can change
        for x = 2:N-1
            for y = 2:N-1
                if change_mask(y,x) == 1
                    % calculate E from surrounding neighbors (p. 422)
                    E = 0;
                    for i = 1:length(neighbors)
                        E = E + A(y+neighbors(i,2),x+neighbors(i,1))*A(y,x);
                    end
                    E = E*(-J);
                    % calculate E from neighbors as is A(y,x) was flipped
                    A_temp = A; A_temp(y,x) = -A_temp(y,x);
                    E_flipped = 0;
                    for i = 1:length(neighbors)
                        E_flipped = E_flipped + A_temp(y+neighbors(i,2),x+neighbors(i,1))*A_temp(y,x);
                    end
                    E_flipped = E_flipped*(-J);
                    % use equation 12.13 to decide if we keep the flip
                    if E > E_flipped
                        P = 1;
                    else
                        P = exp((E-E_flipped)/T);
                    end
                    if P > rand()
                        A_new(y,x) = -A(y,x);
                    end
                end
                % update the total energy of the lattice (all x and y, not just possible flip locations)
                E = 0;
                for i = 1:length(neighbors)
                    E = E + A(y+neighbors(i,2),x+neighbors(i,1))*A(y,x);
                end
                E = E*(-J);
                E_tot = E_tot + E;
            end
        end
        % calculate state magnetism
        Mt(t) = n * sum(sum(A));
        % set lattice for next iteration
        A = A_new;
        % calculate average E
        Et(t) = E_tot / (N-2)^2 / 2; 
        % the /2 in here is from 1, it makes the plot look like the one in the textbook 
        % and 2, I think without it we would be double counting the energy between any two neighbors
    end
end

