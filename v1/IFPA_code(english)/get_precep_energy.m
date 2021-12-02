%Node energy consumption calculation
function [area_radio,radio_area_arr] =  get_precep_energy(sersor_r)
    [~,sersor_size] = size(sersor_r);%Get the number of nodes corresponding to an individual
    radio_area_arr = zeros(sersor_size,1);%Store each individual's perceived energy consumption
    for i=1:sersor_size
        radio_area_arr(i,1) = pi* sersor_r(1,i)^2;
    end
    area_radio = sum(radio_area_arr);%Total radiation area
end