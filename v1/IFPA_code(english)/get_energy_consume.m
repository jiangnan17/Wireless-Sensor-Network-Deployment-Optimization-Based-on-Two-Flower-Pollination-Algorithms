% Energy consumed,formula 5-12 in the paper
% 1. Travelling distance consumption 1mw / m
% 2. Consumption of receiving and sending  mw
% 3. The consumption of radiated area is fixed. 1mw / m2
function [energy,energy_rate] =  get_energy_consume(first_init_wolf,any_indivi,serson_r,energy_init)
    
    energy_ideal = sum(energy_init);%Get the sum of the initial energy
    [~,value] = get_match_value(first_init_wolf,any_indivi,serson_r);
    % Get all distances assigned
    distance = value;%Get the total distance cost
    energy1 = distance * 0.0002;

    

  
   num_send_receive =  get_flood_protocol(any_indivi,serson_r);%
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));
  
    
   [area_radio,~]  = get_precep_energy(serson_r);
   energy3 = area_radio * 0.0003;%Energy consumed by radiation (perception)
   
    
    


    energy = energy1 + energy2 + energy3;
    
    
    energy_rate = energy/energy_ideal;
   
end