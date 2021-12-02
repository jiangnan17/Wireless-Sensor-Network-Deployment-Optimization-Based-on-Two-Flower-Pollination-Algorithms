%Store the matching solutions of storage nodes by LAPJV algorithm
function [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,any_indivi,sensor_r)

    [~,N] = size(first_init_wolf);% 
    
    num_7 = 0;%半径为7的个数
    num_6 = 0;
    num_5 = 0;
    
    %Count the number of three types of nodes
    for i=1:N
        if sensor_r(1,i) == 7
            num_7 = num_7 + 1;
        elseif sensor_r(1,i) == 6
            num_6 = num_6 + 1;
        else
            num_5 = num_5 + 1;
        end
    end
    
    num_7_sensor = (1:1:num_7);%The first type corresponds to the serial number   
    num_6_sensor = ((num_7+1):1:(num_7 + num_6));
    num_5_sensor = ((num_7 + num_6 + 1):1:N);

    matrix_dis_first_any = zeros(N,N);% 
    % Is a good distance, the distance between a node and each node in another individual, if it does not belong to the same category, the distance is set to inf
    for i=1:N
        for j=1:N
            if ((any(num_7_sensor == i))==1) && ((any(num_7_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_7_sensor == i))==1) && ((any(num_7_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf;
            end
            
            if ((any(num_6_sensor == i))==1) && ((any(num_6_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_6_sensor == i))==1) && ((any(num_6_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf; 
            end
            
            
            if ((any(num_5_sensor == i))==1) && ((any(num_5_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_5_sensor == i))==1) && ((any(num_5_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf; 
            end
            
        end
    end

    % save matrix_dis_first_best.mat matrix_dis_first_best;
    match = zeros(N,2);%The first column stores the start, the second column stores the destination
    match(:,1) = (1:1:N);%Solve the order of the first column
    
    %disp(matrix_dis_first_any);
    [match(:,2),value] = lapjv(matrix_dis_first_any);%Call the lapjv algorithm

end
