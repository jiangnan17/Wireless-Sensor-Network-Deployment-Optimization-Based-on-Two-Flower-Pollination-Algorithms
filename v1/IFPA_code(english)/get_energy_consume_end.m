% Total energy consumed when IFPA optimization stops,formula 5-12 in the paper
% 1. Travelling distance consumption 1mw / m
% 2. Consumption of receiving and sending  mw
% 3. The consumption of radiated area is fixed. 1mw / m2
function [match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(first_init_wolf,best_indivi,sersor_r,energy_init)
    
    energy_ideal = sum(energy_init);%Get the sum of the initial energy
    [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,best_indivi,sersor_r);
    % Get all distances assigned
    distance = value;% 
    energy1 = distance * 0.0002;
     
    
    [~,N] = size(first_init_wolf);%
    dis_match_first_best = zeros(N,2);%Save the distance between points after matching
    dis_match_first_best(:,1) = (1:1:N);%The first column is the order, used to store the distance moved by each node after matching
    for i=1:N
        dis_match_first_best(i,2) = matrix_dis_first_any(match(i,1),match(i,2));
    end

    
   % Construct an array like this. The first column is the serial number of the node. The second column is the transmission. The third column is the reception. The fourth column is the fixed noise energy consumption.
   num_send_receive =  get_flood_protocol(best_indivi,sersor_r);%
   
   energy_send_receive_match = zeros(N,2);% Energy consumption of sending and receiving after matching
   energy_send_receive_match(:,1) = num_send_receive(:,1);
   for i=1:N
       energy_send_receive_match(i,2) = 0.0005*sum(num_send_receive(i,2)) + 0.0005*sum(num_send_receive(i,3)) + 0.00005* sum(num_send_receive(i,4));
   end
   
   
   %Do a restore match
   energy_send_receive_match_temp = energy_send_receive_match;
   for i=1:N
       energy_send_receive_match(i,2) = energy_send_receive_match_temp(match(i,2),2);
   end
   
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));

   
   [area_radio,radio_area_arr]  = get_precep_energy(sersor_r);
   energy3 =area_radio * 0.0003;

    
  % Perceived energy consumption after matching, in fact, there is no need to match, because they are all 1
   energy_precep_match = zeros(N,2);
   energy_precep_match(:,1) = (1:1:N);
   energy_precep_match(:,2) = radio_area_arr;
    


    energy = energy1 + energy2 + energy3;
     
    
    energy_rate = energy/energy_ideal;
    
end