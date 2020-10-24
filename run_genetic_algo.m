function [new_gen_2] = run_genetic_algo(gen_count, X1, global_min_x1, global_min_x2, U, L, bit_count, n_feature, N, seed_tab, plot_res_1, plot_res_2, show_tables)
    
    % Step 0: Check if we have more than N samples.
    % Bring input list back to N number of rows
    
    % Count number of rows
    row_in = size(X1, 1);
    if (row_in > N)
        [X1] = bring_back_N(X1, n_feature);
    end
    
    % Step 0.1
    num_to_gen = N * n_feature;
    rand_ls = random_generator(num_to_gen, 0, 0.98); % Generate prob between [0, 0.98]
    
    % Step 1: Evalute objective function and find cumulative fit
    [fit_g1, cumu_fit_g1] = eval_obj(X1);
    
    % Step 2: Find the most fit individual and remember it (Elitisim -- Simply copy and
    % paste)
    
    % TO-DO: Add a user choice here later on
    fit_mat = [X1, transpose(fit_g1)];
    most_fit_g1 = find_most_fit(fit_mat);
    
    % Step 3 Evaluate fraction of fitness and cumulative probability
    [frac_fit_g1, cumu_prob_g1] = eval_fraction_fitness(fit_g1, cumu_fit_g1);
    
    % Separate out x1 and x2 values for all candidates
    % Column vectors here
    x1 = X1(:,1);
    x2 = X1(:,2);
    % Transpose to row vector
    x11 = x1';
    x22 = x2';
    
    % Step 3a Form data table for seed generation
    if(seed_tab == 1)
        candidate_index = transpose([1:N]);
        fit = transpose(fit_g1);
        fractional_fit = transpose(frac_fit_g1);
        cumulative_prob = transpose(cumu_prob_g1);
        if (show_tables == 1)
            % First Table
            DataTable = table(candidate_index, x1, x2, fit,fractional_fit, cumulative_prob)
            % Second table
            fitness_sum = cumu_fit_g1;
            fitness_average = fitness_sum / N;
            ResultTable = table(fitness_sum, fitness_average)
        end
    else
        candidate_index = transpose([1:N]);
        fit = transpose(fit_g1);
        fractional_fit = transpose(frac_fit_g1);
        cumulative_prob = transpose(cumu_prob_g1);
        if (show_tables == 1)
            % First Table
            DataTable = table(candidate_index, x1, x2, fit,fractional_fit, cumulative_prob)
            % Second table
            fitness_sum = cumu_fit_g1;
            fitness_average = fitness_sum / N;
            ResultTable = table(fitness_sum, fitness_average)
        end
    end
        
    % Step 3c Plot top-down location of design variables
    fname = append("Generation",' ',num2str(gen_count));
    if (plot_res_1 == 1)
        figure; % Initiates a separate figure window
        scatter(x11,x22,'k','filled'); % Black star, design candidates
        grid on
        hold on;
        scatter(global_min_x1,global_min_x1,'dr','filled'); % Red diamond indicates target
        xlabel("x1");
        ylabel("x2");
        title(append(fname,' ',"2D plot"));
        xlim([L,U]);
        ylim([L,U]);
        hold off;
        grid off
    end
    
    % Step 3d Plot fitness vs (x1,x2) 3D plot
    if (plot_res_2 == 1)
        z = transpose(fit_g1);
        figure;
        scatter3(x11,x22,z, 'filled', 'MarkerEdgeColor','k');
        view([-37.36 18.64]);
        xlim([L,U]);
        ylim([L,U]);
        zlim([0, 4000]); % HARDCODED zlim
        xlabel("x1");
        ylabel("x2");
        zlabel("fitness");
    end
    
    %----------------------------------------------------------------
    
    % Step 4 Select N random numbers from pool for genetic cross over
    % datasample -- MATLAB's builtin function
    comp_rand = datasample(rand_ls, N,'Replace',false);
    
    % Step 5 Determine mate using wheel of fortune
    mate_matrix = find_mates(X1, cumu_prob_g1, comp_rand);
    
    % Step 6 Reshape to have candiates form couple per row
    mmx = reshape_long_row(mate_matrix);
    
    % Step 7 Perform cross over

    % How many couples do we have
    mmx_couple_count = size(mmx,1);
    
    % Create (N/2) random numbers, row vector
    cross_rand = random_generator(mmx_couple_count, 0, 0.98);
    
    % Perform cross over
    new_gen_2 = create_new_gen(mmx,bit_count,U,L,n_feature);
    
    % Step 8 Reshape new generation matrix to N by n_feature candiates
    new_gen_2 = reshape_N_by_n_feature(new_gen_2, N, n_feature);
    
    % Step 9 Perform elitism, copy-paste the strongest guy from the last generation
    new_gen_2 = [new_gen_2;most_fit_g1];
    
end

% --------------------------------------------------------------------------

function [X1_best] = bring_back_N(X1, n_feature)
    % Calculate fitness and take out the weakest member
    [fit_g1, cumu_fit_g1] = eval_obj(X1);
    fit_g1 = transpose(fit_g1); % Turn row into column vector
    X3 = [X1, fit_g1]; % (N by 3)
    [M,I] = min(X3(:,3));
    X3_1 = X3(1:(I-1),:);
    X3_2 = X3((I+1):end,:);
    
    % Drop (n_feature + 1) column we will add it back later
    X3_1 = X3_1(:,(1:n_feature));
    X3_2 = X3_2(:,(1:n_feature));
    X3 = [X3_1;X3_2];
    
    % Return candidates with N rows
    X1_best = X3;
    
end

% ---------------------------------------------------------------------------

function [fitness, cumu_fit]  = eval_obj(X) % X is a matrix containing x1,x2 . 
    % This function must be modified on case by case basis
    % Unpack and form row vectors for x1 and x2
    x1 = X(:,1);
    x2 = X(:,2);
    
    % Element wise operation - part of vectorized method of evaluating a
    % function with multiple design values instead of using for loop
    % MATLAB's sin, cos automatically uses radian
    obj_f = 100 .* (x2 - x1.^2).^2 + (1 - x1).^2;
    obj_f = transpose(obj_f); % Converting to a row vector
    
    % Remember we have to find fitness by using this equation fit = 1.5 - f
    % where f = obj_f
    
    % The '4000' scalar was found empirifaclly
    fitness = 4000 - obj_f;
    
    % Add up all fitness values
    cumu_fit = sum(fitness,2); % Value 2 tells MATLAB to sum along all the rows
end

% ---------------------------------------------------------------------------

function [most_fit] = find_most_fit(fit_mat) 
    B = sortrows(fit_mat,3, 'descend'); % Most fit individual on top
    most_fit = B(1,[1:2]);
end

% -----------------------------------------------------------------------------

function [frac_to_fit, cumu_prob] = eval_fraction_fitness(fitness, cumu_fit)
    frac_to_fit = fitness./cumu_fit;
    cumu_prob = cumsum(frac_to_fit,2);
end

% ---------------------------------------------------------------------------

function X22 = find_mates(X_current, cumu_prob, comp_rand)
    % Note, both cumu_prob and comp_rand are row vectors
    % The are needed to be flipped to column vector
    
    X11 = [X_current, (cumu_prob).', (comp_rand).']; % Big matrix
    NX = size(X11,1); % Number of rows i.e number of candiates
    X22 = zeros(NX,2);
    num_to_iter = size(X11,1);

    % Scratch pad variable
    row_prob = 0;
    rand_prob = 0;
    test_idx = 1;
    cnt = 1;

    % Counter to keep track of how many counts we did
    for ii = 1: num_to_iter
        row_prob = X11(ii,3); % Choose the cumu_probability for this candidate
        %fprintf("\nCumu prob now --> %f\n",row_prob);
        %fprintf("------------------------\n");
        
        % Look thru each "random" probability sequentially
        for kk = 1: num_to_iter
            rand_prob = X11(kk,4); % Get probability for position at kk
            %fprintf("Rand prob now --> %f\n",rand_prob);
            %%%Test if current cumu probability is greater than current value of random probability
            if (row_prob >= rand_prob)
                X22(kk,:) = X_current(ii,:);
                X11(kk,4) = 1.1; % This will prevent this position random number to evaluated again in next round
                
                %fprintf("-----------------------------------------------------------\n");
                %fprintf("Assigned G_%d position %d going to next candidate\n", ii,kk);
                %fprintf("-----------------------------------------------------------\n");
                break % We stop the inner loop and move onto next candidate
            else
                cnt = cnt + 1; % Dummy to prevent error
                if (cnt<num_to_iter)
                    continue % Go to next random probability value
                else
                    % change actual probability value of this candidate to 1
                    % Exceeded num_iter count
                    % Reset cnt
                    % reset kk
                    % continue again
                    row_prob = 0.99999;
                    cnt = 1;
                    kk = 1;
                    continue
                end
            end
        end
    end
% X22 now is in N x n_feature form. We will reshape them into (N/2,
% n_feature *2) for easier computation of next generation design candidates
end

% ---------------------------------------------------------------------------
function mxx_reshape = reshape_long_row(mate_matrix)
% Function which reshapes N by n_feature to (N/2) by (n_feature*2) matrix
mmx = mate_matrix;
mmx = mmx'; % Shape --> (n_feature by N)
mmx = reshape(mmx,1,[]); % All elements flattened into a row vector
num_count = size(mmx,2);
szz1 = size(mate_matrix,1);
szz2 = size(mate_matrix,2);
ppx = zeros((szz1/2),(szz2 * 2)); % Temporary array to hold (N/2) by (n_feature*2)

% Record number of rows and columns of the reshape matrix
ppx_row = size(ppx,1);
ppx_col = size(ppx,2);
ctn = 1;
  % This nested loop is actually extracting the genes's for each candidate
  % from the flattened matrix
  for pp_row = 1: ppx_row % Row selector
      for pp_col = 1:ppx_col % column selector
          ppx(pp_row, pp_col) = mmx(1,ctn); 
          ctn = ctn + 1;
       end
  end
mxx_reshape = ppx;
end

% ---------------------------------------------------------------------------

function [mmx_reshape] = reshape_N_by_n_feature(new_gen, N, n_feature)
    % Reshape new_gen matrix to 6 by 2
    tempo_mat = zeros(N,n_feature);
    
    % Get number of rows for new_gen matrix
    % Get number of columns for new_gen matrix
    rr = size(new_gen,1);
    cc = size(new_gen,1);
    feature_choose = 1:n_feature;
    row_idx = 0; % Count up to N 
    for i = 1:rr
        row_sel = new_gen(i,:); % One row, all colums
        row_idx = row_idx + 1; 
        
        % Top row, choose col 1 and 2
        for k = feature_choose
            tempo_mat(row_idx,k) = row_sel(1,k);
        end
        
        row_idx = row_idx + 1; 
        % Botttom row, choose col 3 and 4
        for j = feature_choose
            tempo_mat(row_idx,j) = row_sel(1,(j+n_feature));
        end
    end
    mmx_reshape = tempo_mat; % Load up N by n_feature conversion
end

%---------------------------------------------------------------
