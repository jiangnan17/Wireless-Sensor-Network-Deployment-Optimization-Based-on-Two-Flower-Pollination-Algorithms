% According to the formula 33 in the paper, find the tent mapping sequence
function x_m =  get_tent(N)
x = rand(1,1);% Get a random decimal 

x_m = zeros(1,N);% Stores the generated sequence

% figure(1);
for i=1:N 
    if x<=0.6
        x = x/0.6;
        x_m(1,i) = x;
    else
        x = (1-x)/0.4;
        x_m(1,i) = x;
    end
end

% plot(x_m,'.','MarkerSize',20);
% hold on;