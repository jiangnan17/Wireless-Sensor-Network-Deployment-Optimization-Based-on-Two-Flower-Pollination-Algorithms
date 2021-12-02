% The monitoring area is a polygon, the obstacle is a diamond, 
% and IFPA is a single-objective optimization algorithm, 
% and its optimization object is network coverage. 
% Although some other optimization goals are involved in the code, 
% such as energy consumption and node radiation overflow rate, 
% the reader does not need to consider this issue, and does not need to change the code.


%%Main program

clc;
clear ;
close all;



global N;%Number of nodes
global M;%Number of grids
global L;%The length of the monitoring area rectangle
global W;%The width of the monitoring area rectangle
global Grid_cen_x;%The x coordinate of the grid center
global Grid_cen_y;%The y coordinate of the grid center
global Grid_cen_x_and_y;%The x and y coordinates of the center point of the grid
global ger;%Maximum number of iterations of IFPA algorithm


%The following data can be modified as needed, such as the number of nodes, radius, etc.

L = 50;%Unit is m
W = 50;

%Assume that the area of a grid is 1 square meter
M = 2500;%Total number of grids
N = 25;%The number of sensor nodes is 25
%There are three perception radii, 7, 6, 5.
r_max = 7;%The node's perception radius is 7
r_mid = 6;
r_min = 5;

%The energy corresponding to the three nodes
energy_max = 100;%Maximum energy
energy_mid = 90;
energy_min = 80;


% per_sersons_radius_type = [r_max,r_mid,r_min];

%The number of various node types
r_max_num = 1;%Number of nodes with a perception radius of 7
r_mid_num = 2;
r_min_num = N - r_max_num - r_mid_num; 

pos_limit = [0, 50];% Coordinate range of monitoring area, search range of algorithm
ger = 10;% The maximum number of iterations,it can be changed according to your needs
p=0.8;%Conversion factors for global pollination and local pollination
sizepop = 50;%Population size
dimension = 2;% Spatial dimension

%Assembled structure type
struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);
struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);% Temporary population
per_sersons_radius_type = [r_max,r_mid,r_min];






% To construct a monitoring area model, see Figure 10 in the paper, first construct a trapezoidal monitoring area
% Trapezoidal monitoring area is composed of rectangular monitoring area and four triangles.
% Find the data of the 8 vertices of the trapezoid
% The 8 vertices of the trapezoid form four straight line equations, and k, b can be obtained by calculation
% You can use the solve function of matlab to speed up the calculation
% If you need a better understanding, I think you can draw a geometric figure, and then calm down to calculate the mathematical formula.

syms x y;%Define variables
%Upper left corner
k1 = 1;
b1 = 35;
x1_up = solve(k1*x+b1==50,x);%The intersection point above the slash in the upper left corner
y1_down = solve(k1*0+b1==y,y);

%Lower left corner
k2 = -1;
b2 = 15;
y2_up = solve(k2*0+b2==y,y);
x2_down = solve(k2*x+b2==0,x);


%Upper right corner
k3 = -1;
b3 = 85;
x3_up = solve(k3*x+b3==50,x);
y3_down = solve(k3*50+b3==y,y);

%Lower right corner
k4 = 1;
b4 = -35;
y4_up = solve(k4*50+b4==y,y);
x4_down = solve(k4*x+b4==0,x);

%Get the following data
point = zeros(8,2);%Store the location information of these points, from left to right, from top to bottom
point(1,:) = [x1_up,50];
point(2,:) = [0,y1_down];
point(3,:) = [0,y2_up];
point(4,:) = [x2_down,0];
point(5,:) = [x3_up,50];
point(6,:) = [50,y3_down];
point(7,:) = [50,y4_up];
point(8,:) = [x4_down,0];


%The structure of the diamond-shaped obstacle model only needs to calculate the four vertices of the diamond
point_diamond = zeros(2,4);%Four points of the diamond, the book sequence is stored clockwise 
%Find the four points of the new rhombus

syms x y;
%Upper left corner
k5 = 1;
b5 = 10;


point_diamond(1,1) = 25;
point_diamond(2,1) = 35;

%Upper right corner
k6 = -1;
b6 = 60;

point_diamond(1,2) = 35;
point_diamond(2,2) = 25;

%Lower right corner
k7 = 1;
b7 = -10;

point_diamond(1,3) = 25;
point_diamond(2,3) = 15;

%Lower left corner
k8 = 1;
b8 = 40;

point_diamond(1,4) = 15;
point_diamond(2,4) = 25;







load struct_pop_public.mat;%Load public population
struct_pops = struct_pop_public;%Get population data

load struct_pops_second.mat;%Load the second population, named "struct_ pops_ second"

load struct_first_init_public.mat%Load initial deployment data
struct_first_init = struct_first_init_public;%Get the data of an individual (a set of deployment plans)


% initial post-deployment drawing using the first population individual for initial drawing
% A grid is 1m2
% Find grid center coordinates

X_mat = (0:1:50);%x matrix
Y_mat = (0:1:50);%y matrix
Grid_cen_x = zeros(1,L/1);% Grid center x coordinate
Grid_cen_y = zeros(1,W/1);% Grid center y coordinate

% The sum of the data of the two grid points before and after divided by 2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end

%% throw the horizontal and vertical coordinates into a two-dimensional matrix
% Is used to rotate the coordinates. The first row puts the x axis, the second row puts the y axis, and the same column puts one point
% And store the first line near the x coordinate first, then store the second line up
% Grid center coordinates
% Counting from bottom to top, the first line is the first line

Grid_cen_x_and_y = zeros(L,W,2);%Two-dimensional space
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);% 1 represents x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);% Put the y coordinate on the second line
    end
end

x_pos = struct_first_init.per(1,:);%X-coordinate of the first individual
y_pos = struct_first_init.per(2,:);%Y-coordinate of the first individual
sersors_r = struct_first_init.radius;%The first individual, the node group



%plot
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('Initial deployment diagram');
hold on;





% Put the coordinates of randomly deployed points into the matrix
sensor_mat = zeros(2,N);
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

% Stores the number of three nodes
per_sersons_num = struct_pops(1).sersons_num;



%%Calculate the initial coverage of the joint coverage rate and node radiation spillover rate
%However, only focus on coverage in IFPA
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['Initially deployed network coverage：',num2str(cover_rate)]);
disp(['Initially deployed node radiation overflow rate：',num2str(waste_rate)]);



%Calculate connectivity. Connectivity during initialization of the first set of nodes
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('Connected');
else
    disp('Disconnected');
end




% initial population history value is infinitesimal, inf is infinity
best_fitness = -inf;                         % Population history best fitness  
struct_best_indivi = struct_pop_per;                 % Save outstanding individuals
struct_best_indivi_fitness = struct('cover_rate',[],'waste_rate',[],'energy_rate',[],'function_rate',[]);% Adaptation value structure type
struct_best_indivi_fitness_all = repmat(struct_best_indivi_fitness,[1 1]);% Pre-allocated memory




ifa_cover_fitness = zeros(sizepop,1);% Network coverage of a flower (individual)
ifa_waste_fitness = zeros(sizepop,1);%Radiation spillover rate of a flower (individual) node
ifa_energy_fitness = zeros(sizepop,1);%The network energy consumption rate of a flower (individual)
ifa_function_fitness = zeros(sizepop,1);% Combination of three functions(objectives)
weight_rate = [1,0,0];% Sets the weight, we only optimize the network coverage, so set it to 1,


%% Group update
iter = 1;
record_ger = zeros(ger, 1);          % Record the best fitness value in each iteration
record_pop_ave = zeros(ger,1);       % Record the average value of population fitness



aver_fitness_current = inf;%Average fitness of population (current)
aver_fitness_previous = inf;%Average fitness of population (previous)
%Algorithm enters iterative calculation
while iter <= ger

    disp('Current number of iterations:');
    disp(iter);
    % Calculate the fitness value of each individual in the population
    % Find the fitness value of each optimization goal
    for k=1:sizepop
        sensor_mat(1,:) = struct_pops(k).per(1,:);
        sensor_mat(2,:) = struct_pops(k).per(2,:);
        [ifa_cover_fitness(k,1), ifa_waste_fitness(k,1)] = get_Grid_cover_unit_and_rate_waste(sensor_mat,struct_pops(k).radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
        [~,ifa_energy_fitness(k,1)] = get_energy_consume(struct_first_init.per,struct_pops(k).per,struct_pops(k).radius,struct_pops(k).energy_init);
        ifa_function_fitness(k,1) = weight_rate(1,1) * ifa_cover_fitness(k,1) + weight_rate(1,2)*(1-ifa_energy_fitness(k,1)) + weight_rate(1,3) * (1 - ifa_waste_fitness(k,1));
    end
    
   
    
    [ifa_function_fitness_sort,order_index] = sort(ifa_function_fitness);% Sort the fitness values
    disp(ifa_function_fitness_sort);% Is used to temporarily print the adaptation value of the sequence number, which is convenient for viewing the adaptation value
    

   % Update the optimal fitness value and the optimal individual
    if ifa_function_fitness(order_index(sizepop,1),1) > best_fitness
        
       % Use the structure to save the corresponding various adaptation values, which is already in order, so it is the last one
        struct_best_indivi_fitness_all.cover_rate = ifa_cover_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.waste_rate = ifa_waste_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.energy_rate = ifa_energy_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.function_rate = ifa_function_fitness(order_index(sizepop,1),1);
       
        % Save the best individual
        struct_best_indivi = struct_pops(order_index(sizepop,1));
       
        % Save the optimal adaptation value
        best_fitness = ifa_function_fitness(order_index(sizepop,1),1);
    end
    
    record_ger(iter,1) = best_fitness;% Save the best fitness value in this generation
    
    sum_fitness = sum(ifa_function_fitness);
    record_pop_ave(iter,1) = sum_fitness/sizepop;% Calculate the average fitness of the population
    
   % Operate with a temporary population
    struct_pops_new = struct_pops;

    
% Perform special processing to adapt to the WSN coverage problem.
% Here is mixed with the idea of Figure 8 in the paper
    num_swap = 1;%Try as much as possible, you can set 2, or other.
    for j=1:sizepop
        rand_index_swap_best = randperm(N,num_swap);
        % Get some nodes from the best individuals
        for z=1:num_swap
            struct_pops_new(j).per(:,rand_index_swap_best(1,z)) = struct_best_indivi.per(:,rand_index_swap_best(1,z));
        end
    end
    
    
    % Nonlinear convergence factor in Figure 6 of the paper
    b = 1-(1-(((ger-iter)/ger).^2)).^0.5;
    
    %The main process of flower pollination and some improvement methods
    for i=1:sizepop 
        fitness_current = ifa_function_fitness(i,1);%Get the fitness value of the current individual
        %It seems wrong here, but the effect is very good, which is difficult to explain, or local pollination needs a larger proportion
        if rand(1,1) < p               %Local pollination
            index_rand = randperm(sizepop,2);% Randomly select two, or other
            epsilon = b*rand(2,N);% Change step
            
            % Get new intermediate
            struct_temp_indivi_lo = struct_pops_new(i);
            index_N_serson = randperm(N,N);
            
            % All nodes are processed
            for j=1:N
                struct_temp_indivi_lo.per(:,index_N_serson(1,j)) = struct_pops_new(i).per(:,index_N_serson(1,j)) + epsilon(:,j) .*(struct_pops_new(index_rand(1,1)).per(:,index_N_serson(1,j))- struct_pops_new(index_rand(1,2)).per(:,index_N_serson(1,j)));
                
                
               % Simultaneous transboundary processing(The node is not in the monitoring area), 2 means 2 dimensions
                for k=1:2
                     if struct_temp_indivi_lo.per(k,index_N_serson(1,j)) < pos_limit(1,1) || struct_temp_indivi_lo.per(k,index_N_serson(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_lo.per(k,index_N_serson(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                % Handling four triangles
                % Upper left
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;% Is used to mark the question of whether to add or subtract when parallel lines are sought. 1 indicates subtraction
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                % Process the second area
                % Lower left
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;% 2 means add
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2);
                     end
                 end


                % Process the third area
                % Upper right corner
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(5,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(6,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3);
                     end
                 end

                % Process the fourth area
                % Lower right corner
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(8,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4);
                     end
                 end

                % Handle diamond-shaped obstacles
                % Draw four areas
                % Handling the upper left bevel
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,4)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,4)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))% Left hypotenuse
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;% Addition% is used to mark the problem of parallel lines, whether it is addition or subtraction 1 means subtraction
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5);% Let it be on the parallel diagonal line
                     end
                 end

                % Process the second area
                % Upper right
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;% 2 means add
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6);
                     end
                 end


                % Process the third area
                % Lower right corner
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,3)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,3)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7);
                     end
                 end

                % Process the fourth area
                % Lower left corner
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,4)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,3)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,3)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,4))
                     if (k8 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b8) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k8 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b8);
                     end
                 end
            end
            
            
            % Deal with the best individual
            lo_swap_num = 1;%This variable can be changed, for example set to 2, or other
            lo_serson_index = randperm(N,lo_swap_num);
            for k=1:lo_swap_num
                struct_temp_indivi_lo.per(:,lo_serson_index(1,k)) = struct_best_indivi.per(:,lo_serson_index(1,k));
            end
           
            
            % Adaptation value obtained by local pollination
            % Get new fitness value
            [temp_lo_cover_fitness, temp_lo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_lo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius,struct_temp_indivi_lo.energy_init);
            fitness_temp_function_lo = weight_rate(1,1) * temp_lo_cover_fitness + weight_rate(1,2)*(1-temp_lo_energy_fitness) + weight_rate(1,3) * (1 - temp_lo_waste_fitness);

            % Determine whether its fitness value exceeds the previous one. If yes, replace it
            if fitness_temp_function_lo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_lo;% To replace
                
            end
         else% For global pollination
            L1 = Levy2(dimension,iter,ger);% Get step
            index_rand = randperm(N,N);
            
            % get new intermediate
            struct_temp_indivi_glo = struct_pops_new(i);
            for j=1:N
                struct_temp_indivi_glo.per(:,index_rand(1,j)) = struct_pops_new(i).per(:,index_rand(1,j)) + L1' .*(struct_pops_new(i).per(:,index_rand(1,j))- struct_best_indivi.per(:,index_rand(1,j)));%得到临时花朵
            
                % Simultaneous cross-border processing 2 means x, y coordinates
                for k=1:2
                     if struct_temp_indivi_glo.per(k,index_rand(1,j)) < pos_limit(1,1) || struct_temp_indivi_glo.per(k,index_rand(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_glo.per(k,index_rand(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                
                %Data processing for four triangles (trapezoidal monitoring area)
                %This is the same as the obstacle handling in the local pollination above
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;% 
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1);
                     end
                 end


                 %
                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2);
                     end
                 end


                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(5,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(6,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3);
                     end
                 end

                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(8,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4);
                     end
                 end

                %Diamond-shaped obstacle handling
                %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,4)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,4)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5) >= struct_temp_indivi_glo.per(2,index_rand(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6);
                     end
                 end


                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,3)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,3)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7);
                     end
                 end

                 %
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,4)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,3)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,3)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,4))
                     if (k8 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b8) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k8 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b8);
                     end
                 end
            
            end
            
           % Use the best individual for processing
          
            glo_swap_num = 2;%Here can be 2, or other
            glo_serson_index = randperm(N,glo_swap_num);
            for k=1:glo_swap_num
                struct_temp_indivi_glo.per(:,glo_serson_index(1,k)) = struct_best_indivi.per(:,glo_serson_index(1,k));
            end
            
            % Get new fitness value
            [temp_glo_cover_fitness, temp_glo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_glo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius,struct_temp_indivi_glo.energy_init);
            fitness_temp_function_glo = weight_rate(1,1) * temp_glo_cover_fitness + weight_rate(1,2)*(1-temp_glo_energy_fitness) + weight_rate(1,3) * (1 - temp_glo_waste_fitness);

            % Determine whether its fitness value exceeds the previous one. If yes, replace it
            if fitness_temp_function_glo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_glo;
                
            end
           
         end
        
        %Partial cross-operation, which is not reflected in the paper
        %Here is a bit similar to Figure 8 in the paper, greedy cross strategy
        if rand(1,1) < 0.4; %0.4 can be changed according to your needs
            rand_swap_p = b*rand(1,1);%Used to exchange information
            flower_rand_one_va = randperm(sizepop,2);
            
            % Cross with a nearby flower
            for kk=1:2
                if i ~= flower_rand_one_va(1,kk)
                    swap_index = flower_rand_one_va(1,kk);
                    break;
                end
            end
            
            % For cross operation
            
            p_swap = b*rand(2,N);
            struct_flower_current_ac = struct_pops_new(i);
            struct_flower_one_ac = struct_pops_new(swap_index);
            ac_swap_num = 1;%The number of nodes used for switching can be set to 2 or other
            ac_swap_index = randperm(N,ac_swap_num);
            
            
            for z=1:ac_swap_num
               
                struct_flower_current_ac.per(:,ac_swap_index(1,z)) = p_swap(:,ac_swap_index(1,z)) .* struct_pops_new(i).per(:,ac_swap_index(1,z)) + (1 - p_swap(:,ac_swap_index(1,z)) ) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
                struct_flower_one_ac.per(:,ac_swap_index(1,z)) = (1 - p_swap(:,ac_swap_index(1,z)) ).*struct_pops_new(i).per(:,ac_swap_index(1,z)) + p_swap(:,ac_swap_index(1,z)) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
            end
            
           % Handle obstacles
            for j=1:ac_swap_num
                 
                for k=1:2
                     if struct_flower_current_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_current_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_current_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                     
                     if struct_flower_one_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_one_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_one_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;% 
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                %Handling of diamond-shaped obstacles
                %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,3)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,4))
                     if (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8);
                     end
                 end
                 
				 
                % Another variant
                % Handling obstacles
                %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                % Handle diamond-shaped obstacles
                
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,3)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,4))
                     if (k8 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b8) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k8 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b8);
                     end
                 end
                
   
            end
            
            
            
            % Get new fitness value
            ac_orign_fitness = ifa_function_fitness(i,1);
            one_ac_orign_fitness = ifa_function_fitness(swap_index,1);
            
            [temp_ac_cover_fitness, temp_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_current_ac.per,struct_flower_current_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_current_ac.per,struct_flower_current_ac.radius,struct_flower_current_ac.energy_init);
            fitness_temp_function_ac = weight_rate(1,1) * temp_ac_cover_fitness + weight_rate(1,2)*(1-temp_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_ac_waste_fitness);

            %judge whether its fitness value exceeds the previous one, if yes, replace
            if fitness_temp_function_ac > ac_orign_fitness
               
                struct_pops_new(i) = struct_flower_current_ac;%replace
                ifa_function_fitness(i,1) = fitness_temp_function_ac;
            end

            % Another cross
            [temp_one_ac_cover_fitness, temp_one_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_one_ac.per,struct_flower_one_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % Individual (search agent) current fitness
            [~,temp_one_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_one_ac.per,struct_flower_one_ac.radius,struct_flower_one_ac.energy_init);
            fitness_temp_function_one_ac = weight_rate(1,1) * temp_one_ac_cover_fitness + weight_rate(1,2)*(1-temp_one_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_one_ac_waste_fitness);

            % Determine whether its fitness value exceeds the previous one, and if so, replace it
            if fitness_temp_function_one_ac > one_ac_orign_fitness
                
                struct_pops_new(swap_index) = struct_flower_one_ac;
                ifa_function_fitness(swap_index,1) = fitness_temp_function_one_ac;
            end
            
        end
    end
   
   
   %The tent map is used to enrich the population. The second population will be used here 
   aver_fitness_previous = aver_fitness_current;
   aver_fitness_current = mean(ifa_cover_fitness);
   %Lack of population diversity,θ=0.0003,θ can be changed according to your needs
   if abs(aver_fitness_current - aver_fitness_previous) < 0.0003
       repalce_node_index = randperm(N,1);% 1 here can be set to 2 or 3, etc.
       replace_aget_index = randperm(sizepop,1);% 1 here can be set to 2 or 3, etc.
       struct_pops_new(replace_aget_index(1,1)).per(:,repalce_node_index(1,1)) = struct_pops_second(replace_aget_index(1,1)).per(:,repalce_node_index(1,1));
   end
   
   struct_pops = struct_pops_new;% Update population
   iter = iter+1;
end

%Write to xlsx table

disp('Data writing ...');
xlswrite('ifa_cover_25serson.xlsx',record_ger,1);



figure(2);
plot(record_ger);
title('Optimal fitness value in each generation')

figure(3);
plot(record_pop_ave);
title('The average fitness value of each population')

figure(4);
draw_circle(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius);
title('Optimized deployment diagram');


figure(5);
draw_match(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius);
title('Assignment diagram after deployment');
hold on;





% Judge network connectivity
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(struct_best_indivi.per,struct_best_indivi.radius);

if is_connec == 1
    disp('Connected');
else
    disp('Disconnected');
end

% Draw network minimum spanning tree
figure(6);
draw_MST(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('Minimum spanning tree');
hold on;

% Bring back important news
[match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius,struct_first_init.energy_init);

% The final position after the last match
end_ifa_match = struct_first_init;

% Match, restore
for i=1:N
    % Location update
    end_ifa_match.per(:,i) = struct_best_indivi.per(:,match(i,2));
    
    end_ifa_match.energy_end(1,i) = end_ifa_match.energy_init(1,i) - (0.0002*dis_match_first_best(i,2) + energy_send_receive_match(i,2) + 0.0003*energy_precep_match(i,2)); 
end
disp(['Total energy during network initialization：',num2str(sum(end_ifa_match.energy_init(1,:)))]);
disp(['Total distance moved after node deployment：',num2str(sum(dis_match_first_best(:,2)))]);
disp(['Energy consumed by moving distance after node deployment：',num2str(0.0002*sum(dis_match_first_best(:,2)))]);
disp(['Total energy consumed after node deployment：',num2str(sum(end_ifa_match.energy_init(1,:)) - sum(end_ifa_match.energy_end(1,:)))]);


figure(7);
draw_circle(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius);
title('Deploy the diagram after the node moves to its location');
hold on;


% Is only for cooperation to obtain the minimum spanning tree
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(end_ifa_match.per,end_ifa_match.radius);

figure(8);
draw_MST(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('The minimum spanning tree at the time of last deployment');
hold on;


disp(['Algorithm-optimized function：',num2str(struct_best_indivi_fitness_all.function_rate)]);
disp(['Corresponding network coverage：',num2str(struct_best_indivi_fitness_all.cover_rate)]);
disp(['Corresponding network energy consumption rate：',num2str(struct_best_indivi_fitness_all.energy_rate)]);
disp(['Corresponding node radiation waste rate：',num2str(struct_best_indivi_fitness_all.waste_rate)]);

% Saved final deployment result
end_ifa_match1 = end_ifa_match;
save end_ifa_match1.mat  end_ifa_match1;
