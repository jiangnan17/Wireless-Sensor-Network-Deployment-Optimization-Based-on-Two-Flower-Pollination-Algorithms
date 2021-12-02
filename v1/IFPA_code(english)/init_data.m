   
%Initialize population data, two populations are generated,
%one of which is used to initialize the population and enrich
%the diversity of the population,use public population data

clc;
clear ;
close all;
format long;



%Delete previously initialized populations and data

delete('struct_pop_public.mat');%Whole population
delete('struct_first_init_public.mat');%An individual in a population


global N;%Number of nodes
global M;%Number of grids
global L;%The length of the monitoring area rectangle
global W;%The width of the monitoring area rectangle
global Grid_cen_x;%The x coordinate of the grid center
global Grid_cen_y;%The y coordinate of the grid center
global Grid_cen_x_and_y;%The x and y coordinates of the center point of the grid

%The following data can be modified as needed, such as the number of nodes, radius, etc.

L = 50;%Unit is m
W = 50;

%Assume that the area of a grid is 1 square meter
M = 2500;%Total number of grids

%There are three perception radii, 7, 6, 5.
r_max = 7;%The node's perception radius is 7
r_mid = 6;
r_min = 5;

%The energy corresponding to the three nodes
energy_max = 100;%最大的能量
energy_mid = 90;
energy_min = 80;



per_sersons_radius_type = [r_max,r_mid,r_min];


N = 25;%The number of sensor nodes is 25
sizepop = 50;%Population size
dimension = 2;% Spatial dimension

%The number of various node types
r_max_num = 1;%Number of nodes with a perception radius of 7
r_mid_num = 2;
r_min_num = N - r_max_num - r_mid_num; 

%Assembled structure type
struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);
struct_pops = repmat(struct_pop_per,[1 sizepop]);% Generate structure array, 1 row and 50 columns


% Is used to store energy, but in IFPA, energy consumption is not considered
energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);

% Initialize the energy and radius of each node This is uniform, 
% although energy consumption in IFPA is not considered.
for k=1:N
    if k>=1 && k <=r_max_num 
       radius_arr(1,k) = r_max; 
       energy_init_arr(1,k) = energy_max;
       energy_end_arr(1,k) = energy_max;
    elseif k>=r_max_num+1 && k <= (r_max_num + r_mid_num)
       radius_arr(1,k) = r_mid; 
       energy_init_arr(1,k) = energy_mid;
       energy_end_arr(1,k) = energy_mid;
    else
       radius_arr(1,k) = r_min; 
       energy_init_arr(1,k) = energy_min;
       energy_end_arr(1,k) = energy_min;
    end
end

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

% Generate nodal data, energy, radius structure
% Simultaneously handle the problem if the node falls into the obstacle
for i=1:sizepop
    
   % Initialize the radius and energy. Because the use of radius is involved below, initialize the radius first
    struct_pops(i).radius(1,:) = radius_arr;
    struct_pops(i).energy_init(1,:) = energy_init_arr;
    struct_pops(i).energy_end(1,:) = energy_end_arr;
    
   % Number of initialized nodes
    struct_pops(i).sersons_num = [r_max_num,r_mid_num,r_min_num];
    
    
    % Initialize node position coordinates,the formula 33 in the paper,
    %another population will be generated at the bottom of this document
    x = L*get_tent(N);%Randomly generate the x coordinate of the node
    y = W*get_tent(N);%Randomly generate the y coordinate of the node
    
    %If the tent mapping effect is not very good, you can use the following method instead
%     x = L*rand(1,N);
%     y = W*rand(1,N);
    
    struct_pops(i).per(1,:) = x;
    struct_pops(i).per(2,:) = y;
    
    %Refer to Figure 5 in the paper for the handling when the% node falls into an obstacle.
    for j=1:N
          % Draw four areas
          % Handle the triangle on the upper left bevel
         if (struct_pops(i).per(1,j)>=0&&struct_pops(i).per(1,j)<=point(1,1)) && (struct_pops(i).per(2,j)>=point(2,2)&&struct_pops(i).per(2,j)<=50)
             if (k1 * struct_pops(i).per(1,j) + b1) <= struct_pops(i).per(2,j)%Left hypotenuse
                 k1 = 1;
                 b1 = 35;
                 case_b = 1;% It is used to mark the question of whether to add or subtract when parallel lines are sought. 1 indicates subtraction
                 [k1,b1] = get_new_function(k1,b1,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k1 * struct_pops(i).per(1,j) + b1);% Let it be on the parallel diagonal line
             end
         end
         
         
         % Process the second area
         % Lower left
         if (struct_pops(i).per(1,j)>=0&&struct_pops(i).per(1,j)<=point(1,1)) && (struct_pops(i).per(2,j)>=0&&struct_pops(i).per(2,j)<=point(3,2))
             if (k2 * struct_pops(i).per(1,j) + b2) >= struct_pops(i).per(2,j)
                 k2 = -1;
                 b2 = 15;
                 case_b = 2;% 2 means add
                 [k2,b2] = get_new_function(k2,b2,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k2 * struct_pops(i).per(1,j) + b2);
             end
         end
         
         
          % Process the third area
          % Upper right corner
         if (struct_pops(i).per(1,j)>=point(5,1)&&struct_pops(i).per(1,j)<=50) && (struct_pops(i).per(2,j)>=point(6,2)&&struct_pops(i).per(2,j)<=50)
             if (k3 * struct_pops(i).per(1,j) + b3) <= struct_pops(i).per(2,j)
                 k3 = -1;
                 b3 = 85;
                 case_b = 1;
                 [k3,b3] = get_new_function(k3,b3,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k3 * struct_pops(i).per(1,j) + b3);
             end
         end
         
        % Process the fourth area
        % Lower right corner
         if (struct_pops(i).per(1,j)>=point(8,1)&&struct_pops(i).per(1,j)<=50) && (struct_pops(i).per(2,j)>=0&&struct_pops(i).per(2,j)<=point(7,2))
             if (k4 * struct_pops(i).per(1,j) + b4) >= struct_pops(i).per(2,j)
                 k4 = 1;
                 b4 = -35;
                 case_b = 2;
                 [k4,b4] = get_new_function(k4,b4,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k4 * struct_pops(i).per(1,j) + b4);
             end
         end
         
         
        % Dealing with the problem of nodes entering diamond-shaped obstacles
        % Draw four areas
        % Handling the upper left bevel
         if (struct_pops(i).per(1,j)>=point_diamond(1,4)&&struct_pops(i).per(1,j)<=point_diamond(1,1)) && (struct_pops(i).per(2,j)>=point_diamond(2,4)&&struct_pops(i).per(2,j)<=point_diamond(2,1))
             if (k5 * struct_pops(i).per(1,j) + b5) >= struct_pops(i).per(2,j)% Left hypotenuse
                 k5 = 1;
                 b5 = 10;
                 case_b = 2;% Addition% is used to mark the problem of parallel lines, whether it is addition or subtraction 1 means subtraction
                 [k5,b5] = get_new_function(k5,b5,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k5 * struct_pops(i).per(1,j) + b5);% Let it be on the parallel diagonal line
             end
         end
         
         
        % Process the second area
        % Upper right
         if (struct_pops(i).per(1,j)>=point_diamond(1,1)&&struct_pops(i).per(1,j)<=point_diamond(1,2)) && (struct_pops(i).per(2,j)>=point_diamond(2,2)&&struct_pops(i).per(2,j)<=point_diamond(2,1))
             if (k6 * struct_pops(i).per(1,j) + b6) >= struct_pops(i).per(2,j)
                 k6 = -1;
                 b6 = 60;
                 case_b = 2;% 2 means add
                 [k6,b6] = get_new_function(k6,b6,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k6 * struct_pops(i).per(1,j) + b6);
             end
         end
         
         
        % Process the third area
        % Lower right corner
         if (struct_pops(i).per(1,j)>=point_diamond(1,3)&&struct_pops(i).per(1,j)<=point_diamond(1,2)) && (struct_pops(i).per(2,j)>=point_diamond(2,3)&&struct_pops(i).per(2,j)<=point_diamond(2,2))
             if (k7 * struct_pops(i).per(1,j) + b7) <= struct_pops(i).per(2,j)
                 k7 = 1;
                 b7 = -10;
                 case_b = 1;
                 [k7,b7] = get_new_function(k7,b7,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k7 * struct_pops(i).per(1,j) + b7);
             end
         end
         
        % Process the fourth area
        % Lower left corner
         if (struct_pops(i).per(1,j)>=point_diamond(1,4)&&struct_pops(i).per(1,j)<=point_diamond(1,3)) && (struct_pops(i).per(2,j)>=point_diamond(2,3)&&struct_pops(i).per(2,j)<=point_diamond(2,4))
             if (k8 * struct_pops(i).per(1,j) + b8) <= struct_pops(i).per(2,j)
                 k8 = 1;
                 b8 = 40;
                 case_b = 1;
                 [k8,b8] = get_new_function(k8,b8,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k8 * struct_pops(i).per(1,j) + b8);
             end
         end
         
    end
    
end



%% initial post-deployment drawing using the first population individual for initial drawing
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


x_pos = struct_pops(1).per(1,:);%X-coordinate of the first individual
y_pos = struct_pops(1).per(2,:);%Y-coordinate of the first individual
sersors_r = struct_pops(1).radius;%The first flower individual, the node group




% plot
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('Initial deployment diagram');
hold on;
struct_first_init_indi = struct_pops(1);%Save the first individual data

% Save a copy of the data of the first individual for sharing
struct_first_init_public = struct_first_init_indi;
save struct_first_init_public.mat struct_first_init_public;



% Put the coordinates of randomly deployed points into the matrix
sensor_mat = zeros(2,N);
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

% Stores the number of three nodes
per_sersons_num = struct_pops(1).sersons_num;


% Initialized to get the joint probability
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['Coverage of initial network deployment：',num2str(cover_rate)]);




%Calculate connectivity. Connectivity during initialization of the first set of nodes
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('Connected');
else
    disp('Disconnected');
end



% Because the first individual in the first generation is used as the initial, it is necessary to consider regenerating the first individual and only generate coordinates
% In order to reduce the amount of calculation, a new individual is pieced together from sizepop individuals

new_first_per = zeros(2,N);
piece_together_index = randperm(N,N);
rand_indivi = randperm(sizepop,N);% Note that N is less than sizepop
for i=1:N
    new_first_per(:,i) = struct_pops(rand_indivi(1,i)).per(:,piece_together_index(1,i));
end
struct_pops(1).per = new_first_per;
struct_pop_public = struct_pops;% Public initial population


%Generate a spare(another) population to make the algorithm jump out of the local optimum
struct_pops_second = struct_pops;
%Only coordinates need to be changed

%Replace some nodes of some search agents in a cross way
replace_aget_index = randperm(sizepop,sizepop);
node_temp = zeros(2,1);%Temporary variable, used to exchange position coordinates of nodes
for i1=1:sizepop/2
    repalce_node_index = randperm(N,2);% 2 here can be set to 4 or 6, etc.
    node_temp = struct_pops_second(replace_aget_index(1,i1)).per(:,repalce_node_index(1,1));
    struct_pops_second(replace_aget_index(1,i1)).per(:,repalce_node_index(1,1)) = struct_pops_second(replace_aget_index(1,sizepop - i1)).per(:,repalce_node_index(1,2));
    struct_pops_second(replace_aget_index(1,sizepop - i1)).per(:,repalce_node_index(1,2)) = node_temp;
end


% Save the population to the wolf_pop_public.mat file
save struct_pop_public.mat struct_pop_public;
save struct_pops_second.mat struct_pops_second;
disp('Data initialization has been completed');


