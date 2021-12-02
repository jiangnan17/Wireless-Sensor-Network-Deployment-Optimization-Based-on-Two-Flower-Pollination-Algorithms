%%Determine whether the network is connected
%According to formula 19 in the paper
function [is_connec,adjacencyMatrix,adjacencyMatrix_dis ]= get_connection(any_indivi,sersors_r)

[~,serson_size] = size(any_indivi);%Get the size of the sensor node


%%Computing connectivity
adjacencyMatrix=zeros(serson_size,serson_size);% Defines the sensor interconnection adjacency matrix, which is a directed graph
adjacencyMatrix_dis=ones(serson_size,serson_size).*inf;
for i=1:1:serson_size
    for j=(i+1):1:serson_size   %Because there are multiple radii, the judgment of network connectivity has to be handled slightly

            

            flag1 = 0;%Indicates that A-> B is not connected
            flag2 = 0;%Indicates that B-> A is not connected
            dis_squ = (any_indivi(1,i)-any_indivi(1,j))^2 +...
        (any_indivi(2,i)-any_indivi(2,j))^2;
           %Add distance
            if dis_squ <= (2 * sersors_r(1,i))^2;%Perceivable between two points
                flag1 = 1;%A-> B one way arrival
            end
            if dis_squ <= (2 * sersors_r(1,j))^2;%Perceivable between two points
                flag2 = 1;%B-> A one way arrival
            end
            
            %Construct an undirected graph, when A-> B B-> A is regarded as an undirected graph.
            if flag1==1 && flag2==1;%；You can feel each other between two points
                adjacencyMatrix(i,j)=1;%
                adjacencyMatrix(j,i)=1;
                
                %Join distance
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end   
    end
end 

%Save the last directed graph adjacency matrix
% save adjacencyMatrix_dis.mat adjacencyMatrix_dis;
% save adjacencyMatrix.mat adjacencyMatrix;

S=zeros(serson_size,serson_size);
for m=1:1:serson_size-1
    S=S+adjacencyMatrix^m;  %全部加到S
end;

%S=M+M^2+M^3+...+M^(N-1)；

is_connec = 0;%Initialize as disconnected
if all(all(S))==1 %Determine whether the S vector is all non-zero elements
    is_connec = 1;%Connect if all are non-zero
end
end