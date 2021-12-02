x = 0.1999;
N_num = 2000;
x_m = zeros(1,N_num);
y_m = zeros(1,N_num);
figure(1);
for i=1:N_num 
    if x<=0.6
        x = x/0.6;
        x_m(1,i) = x;
    else
        x = (1-x)/0.4;
        x_m(1,i) = x;
    end
%     y_m(1,i) = i*(1/N_num);
    y_m(1,i) = i;
end
plot(y_m,x_m,'.','MarkerSize',20);
hold on;