function new_gen = create_new_gen(mmx,bit_count,U,L,n_feature)
    % Initialize working variables
    row_idx = size(mmx,1);
    col_idx = size(mmx,2);
    cc_mat = zeros(row_idx, col_idx); % Temporary matrix to hold children in 3 by 4 shape
    ccx1 = [];
    ccx2 = [];
    % Create (N/2) random numbers
    cross_rand = random_generator(row_idx, 0, 0.98); % (N/2) number of cross over points
    
    % Select one row, fill in columns sequentially
    for i = 1:row_idx
        %fprintf("Cross over iteration %d\n", i);
        in_mmx = mmx(i,:);
        in_cross_rand = cross_rand(:,i);
        [ccx1,cxx2] = genetic_crossover(in_mmx, in_cross_rand,bit_count,U,L,n_feature);
        cons_mat = [ccx1, cxx2]; % Convert two row vector in 1
        % Create new_gen matrox
        for k = 1: col_idx
            cc_mat(i,k) = cons_mat(1,k);
        end
        new_gen = cc_mat; % Shape (3 by 4) for N =6 and n = 2
    end
    % Now we reshape 3 by 4 matrix into 6 by 2 to match out table
    
end

% ------------------------------------------------
function [child1, child2] = genetic_crossover(couple_row, rand_num, bit_count,U,L,n_feature)
    % Length of one 'chromosome' in bits
    nn_chromo = bit_count * n_feature; % if 8 bits and 2 genes then one candidate -- 16 bits 
    J = 2^bit_count - 1;
    
    % Define ranges for bit level manipulation
    bit_range = create_bit_range_mat(bit_count, n_feature);
    
    % --------------- HARDCODED SECTION FOR TWO VARIABLE ------------- %
    
    % Encode scheme
    % variable -> base_10 -> bit_count binary -> nn_chromo bit parent
    
    % Unpack row vector into four gene sequences
    mx1 = couple_row(:,1);
    mx2 = couple_row(:,2);
    fx1 = couple_row(:,3);
    fx2 = couple_row(:,4);
    
    % Convert variable to X10
    mx1 = VartoX10(mx1,U,L,J);
    mx2 = VartoX10(mx2,U,L,J);
    fx1 = VartoX10(fx1,U,L,J);
    fx2 = VartoX10(fx2,U,L,J);
    
    % Cross-over point calculation needs to be updated
    % Generate one cross over point
    N= nn_chromo - 1; 
    CrossoverIndex = cast(rand_num * N,'uint16');
    %fprintf("Cross-over point %d\n\n", CrossoverIndex);
    
    % Encode decimal to binary
    mx1 = de2bi(mx1, bit_count);
    mx2 = de2bi(mx2, bit_count);
    fx1 = de2bi(fx1, bit_count);
    fx2 = de2bi(fx2, bit_count);
    
    % Make nn_chromo bit genome sequence (each feature placed side by side)
    male_gen = [mx1, mx2];
    female_gen = [fx1, fx2];
    
    % Print out 16 bit male and female gene sequence
    %disp("Male"); print_genome(male_gen)
    %disp("Female"); print_genome(female_gen)
    
    % Crossover, binary values
    child1 = [male_gen(1:CrossoverIndex) female_gen(CrossoverIndex+1:end)];
    child2 = [female_gen(1:CrossoverIndex) male_gen(CrossoverIndex+1:end)];
    
    %---------------------------- Mutation operator code -----------------
    
    % Bit-Flip Mutation
    % Define ranges
    min_pm = 0.001; % Usual minimum mutation probability
    %max_pm = 1/nn_chromo; % Ref - 3
    max_pm = 0.006; % Ref - 3
    
    % Generate nn_chromo mutation probabilities between min_pm to max_pm
    prob_mutation = random_generator(nn_chromo,min_pm,max_pm);
    % Shuffle to increase chance of randomization
    prob_mutation = randomize_array(prob_mutation);
    
    % Pick one probability randomly
    rng shuffle; % Reset random seed generator
    prob_mutation = randsample(prob_mutation, 1, 'true', prob_mutation);
    
    % Create a random probability for each bit location
    bit_mutation_prob = random_generator(nn_chromo,min_pm,max_pm);
    
    % Shuffle to increase chance of randomization
    bit_mutation_prob = randomize_array(bit_mutation_prob);
    
    % Choose a random bit location
    rng shuffle;
    bit_loc = linspace(1,nn_chromo,nn_chromo);
    bit_loc = randsample(bit_loc,1);
    bit_prob = bit_mutation_prob(1,bit_loc);
    
    % HARDCODED for two variable problem
    % Is bit_prob greater than prob_mutation? Yes, flip that bit
    if (bit_prob > prob_mutation)
        child1(1,bit_loc) = ~child1(1,bit_loc);
        child2(1,bit_loc) = ~child2(1,bit_loc);
    else
        %fprintf("No mutation!\n");
    end 
    
    % --------------------------------------------------------------------
    % Two point SWAP Mutation operator (for later)
    % Code goes here
    
    %---------------------------- Mutation operator code -----------------
    
    % Decode scheme
    % variable <- base_10 <- bit_count binary <- nn_chromo bit parent
    % MATLAB's bit get only works on an integer
    [cx1, cx2, dx1, dx2] = X32toBin(child1, child2, bit_range); 
    
    cx1 = BintoVar(cx1, U,L,J);
    cx2 = BintoVar(cx2, U,L,J);
    dx1 = BintoVar(dx1, U,L,J);
    dx2 = BintoVar(dx2, U,L,J);
    
    % --------------- HARDCODED SECTION FOR TWO VARIABLE ------------- %
    
    % Create the 'children' row vectors
    child1 = [cx1, cx2]; 
    child2 = [dx1, dx2];
end

% Dependency for genetic_crossover

% ------------------------------------------------
function [out_vector] = randomize_array(in_vector)
% Function to randomize location of elements in a row vector
% Test if in_vector is a row vector
V = isrow(in_vector);
    if (V == 0)
        in_vector = transpose(in_vector); % Make sure the in vector is a row vector
    end
rand_pos = randperm(length(in_vector)); %array of random positions
% new array with original data randomly distributed 
    for k = 1:length(in_vector)
        out_vector(k) = in_vector(rand_pos(k));
    end
% If our priginal vector was a column vector, we need to put it back in
% correct shape
    if (V == 0)
        out_vector = transpose(out_vector); % Revert it back
    end
end

% ------------------------------------------------
% Xreal -- value in scale defined for the design variable
function [Xreal] = BintoVar(Xbin, U,L,J)
    Xreal = bi2de(Xbin); 
    % For some reason this has to be explicitly casted to double
    Xreal = double(Xreal);
    % Scale base-10 value to design variable scale
    Xreal = X10toVar(Xreal,U,L,J); 
end


% ------------------------------------------------
function [Xreal] = X10toVar(X10,UU,LL,JJ)
  Xreal = LL + ((UU - LL)/JJ) * X10;
  % MATLAB cannot convert a float into binary
  Xreal = double(Xreal);
end

% ------------------------------------------------
% X10 must be fininte, non-negative integer
function [X10] = VartoX10(XVar,U,L,J)
  X10 = ((XVar - L) * J)/(U - L);
  X10 = int16(X10);
end

%-----------------------------------------------

function [cx1, cx2, dx1, dx2] = X32toBin(child1, child2,bit_range)
    % HARDCODED, needs to change
    % Function to automate MATLAB's way of finding and extracting binary bits
    % from an integer number
    child1 = bi2de(child1);
    child2 = bi2de(child2);
    first_8 = bit_range(1,:); % First row, all columns
    second_8 = bit_range(2,:); % 2nd row, all columns
    cx1 = bitget(child1, first_8);
    cx2 = bitget(child1, second_8);
    dx1 = bitget(child2, first_8);
    dx2 = bitget(child2, second_8);
end

% ----------------------------------------------
function print_genome(P)
    fprintf('Genome : [');
    fprintf('%g ', P);
    fprintf(']\n');
    P = bi2de(P);
    %fprintf("Genome in decimal --> %d\n\n", P);
end
% -----------------------------------------------
function [bit_range_matrix] = create_bit_range_mat(bit_count, n_feature)
    % Function which creates a matrix containing integer locations for
    % bitwise manipulation
    nn_chromo = bit_count * n_feature;
    mmx = linspace(1,nn_chromo,nn_chromo);
    ppx = zeros(n_feature, bit_count); % Temporary matrix
    ctn = 1; % Variable to loop through each element in a row matrix
    % Use the nested loop trick to reshape
    for pp_row = 1: n_feature % Row selector
        for pp_col = 1:bit_count % column selector
            ppx(pp_row, pp_col) = mmx(1,ctn); 
            ctn = ctn + 1;
        end
    end
   bit_range_matrix = ppx; % Output
end
