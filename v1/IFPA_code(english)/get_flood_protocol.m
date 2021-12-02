%Process of simulating flooding protocol
function num_send_receive =  get_flood_protocol(any_indivi,sersor_r)
 [~,sersor_size] = size(any_indivi);

 num_send_receive = zeros(sersor_size,4);
 num_send_receive(:,1) = (1:1:sersor_size);% 
 num_send_receive(:,2) = 1;% 
 num_send_receive(:,4) = 1;% 
 connect_num = zeros(sersor_size,2);%Two columns, the first column corresponds to the serial number of the node, the second column corresponds to the number of neighboring nodes
 connect_num(:,1) = (1:1:sersor_size);%
 %What we receive is equivalent to finding the number of its neighbors, and what we receive is equivalent to what its neighbors can send.
 for i=1:sersor_size
    for j=1:sersor_size%Find distance from nodes other than itself
        if(i==j)
            continue;
        end
        %There should be improvements here a-> b b-> a The communication distance is different whether it is covered or not, now unified as ri + rj
        if ((any_indivi(1,i) - any_indivi(1,j))^2 + (any_indivi(2,i) - any_indivi(2,j))^2)...
                ^0.5 <= 2*sersor_r(1,j);
            connect_num(i,2) = connect_num(i,2) + 1;
        end
    end
 end
    
    num_send_receive(:,3) = connect_num(:,2); %Get the number of neighbor nodes
end