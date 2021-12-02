%%主程序
clc;
clear ;
close all;
%删除相应的文件

global N;
global M;
global L;
global W;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global ger;

p=0.8;%判断是否是全局优化还是局部优化

L = 50;%长
W = 50;%宽
%假设1平方米一个网格
M = 2500;%网格总数
r_max = 7;%感知半径为5
r_mid = 6;
r_min = 5;
energy_max = 100;%最大的能量
energy_mid = 90;
energy_min = 80;

per_sersons_radius_type = [r_max,r_mid,r_min];
%假设大、中为5，剩下为小
N = 25;%30个传感器节点
sizepop = 200;%种群规模
dimension = 2;% 空间维数  前行放x、y，第三行放半径
ger = 2;% 最大迭代次数
pos_limit = [0, 50];            % 设置位置参数限制
%个数限制
r_max_num = 1;%序号为1-5
r_mid_num = 2;%序号为6-10
r_min_num = N - r_max_num - r_mid_num;%序号为11-N

aver_fitness_ger = zeros(3,ger);%存三种适应值每代的平均值
struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);%结构体类型

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%临时的一个种群

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);


%求出梯形的四个点
syms x y;%先定义一个变量
%左上角
k1 = 1;
b1 = 35;
x1_up = solve(k1*x+b1==50,x);%左上角的斜线的上个交点
y1_down = solve(k1*0+b1==y,y);

%左下角
k2 = -1;
b2 = 15;
y2_up = solve(k2*0+b2==y,y);
x2_down = solve(k2*x+b2==0,x);


%右上角
k3 = -1;
b3 = 85;
x3_up = solve(k3*x+b3==50,x);
y3_down = solve(k3*50+b3==y,y);

%右下角
k4 = 1;
b4 = -35;
y4_up = solve(k4*50+b4==y,y);
x4_down = solve(k4*x+b4==0,x);

%以下数据验证完毕，完全正确
point = zeros(8,2);%存储这些点  从左  从上往下
point(1,:) = [x1_up,50];
point(2,:) = [0,y1_down];
point(3,:) = [0,y2_up];
point(4,:) = [x2_down,0];
point(5,:) = [x3_up,50];
point(6,:) = [50,y3_down];
point(7,:) = [50,y4_up];
point(8,:) = [x4_down,0];

load struct_pop_public.mat;%加载该种群
struct_pops = struct_pop_public;%得到种群数据
struct_pops_new_origi = repmat(struct_pop_per,[1 2*sizepop]);%包含原有种群和新的种群

struct_pops_archive = struct_pops_new_origi;%用于保存archive里面的个体


struct_pops_archive_temp = struct_pops_new_origi;%用于保存archive_temp里面的个体


load struct_first_init_public.mat%加载最开始的一个个体数据
struct_first_init = struct_first_init_public;%得到初始化个体数据


%%初始的部署后画图  拿第一个粒子拿去初始画图

%求网格中心坐标
X_mat = (0:1:50);%x矩阵
Y_mat = (0:1:50);%y矩阵
Grid_cen_x = zeros(1,L/1);%网格中心点x坐标
Grid_cen_y = zeros(1,W/1);%网格中心点y坐标
%前后两者相加之和除以2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end


%%把横纵坐标丢到一个二维矩阵当中
%用于转坐标  第一行放x轴 第二行放y轴，同一列放一个点坐标
%且先存靠近x坐标的第一行，然后往上存第二行
%网格中心坐标
Grid_cen_x_and_y = zeros(L,W,2);%共2500个网格中心，但是每个网格中心有x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1代表x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%把y坐标放到第二行
    end
end


x_pos = struct_first_init.per(1,:);%第一个粒子   粒子即是解  的x坐标
y_pos = struct_first_init.per(2,:);%第一个粒子   粒子即是解  的y坐标
sersors_r = struct_first_init.radius;%第一个粒子 




%进行画图
% figure(1);
% draw_circle(x_pos,y_pos,sersors_r);
% title('初始化部署图');
% hold on;





%把随机部署的点的坐标放到矩阵当中
sensor_mat = zeros(2,N);%预分配内存
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

%存放三种节点的数目
per_sersons_num = struct_pops(1).sersons_num;
%存放三种节点的半径


%%初始化得到联合概率和节点浪费率
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['初始化的覆盖率：',num2str(cover_rate)]);
disp(['初始化的浪费率：',num2str(waste_rate)]);



%%计算连通性  第一个节点初始化时的连通性
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('连通');
else
    disp('非连通');
end




%%初始化种群历史值为  无穷小   inf为无穷大

best_fitness = -inf;                         % 种群历史最佳适应度  
struct_best_indivi = struct_pop_per;                 % 保存优秀个体
struct_best_indivi_fitness = struct('cover_rate',[],'waste_rate',[],'energy_rate',[],'function_rate',[]);%适应值结构体类型
struct_best_indivi_fitness_all = repmat(struct_best_indivi_fitness,[1 1]);%预分配内存
%%在上面已经画图了
%%以上为画出初始化时，50个粒子开始的位置


archive = zeros(1+6,2*sizepop);%和pareto设置是一样的  4行是拥挤距离
archive_num = 0;%真正实在的非支配
archive_num_limit = sizepop;%最多领导者的限制
%由于无穷大无法计算  设置一个很大的数字 inf_value  解决轮盘赌法那里的问题


origi_new_ifa_cover_fitness = zeros(2*sizepop,1);%个体覆盖率
origi_new_ifa_waste_fitness = zeros(2*sizepop,1);%个体浪费率
origi_new_ifa_energy_fitness = zeros(2*sizepop,1);%个体能耗率  化成1/xx形式 

ifa_function_fitness = zeros(2*sizepop,1);%三个函数的结合
ifa_fitness_temp = zeros(2*sizepop,1);%临时保存当前适应度
%% 群体更新
iter = 1;
record_ger = zeros(ger, 1);          % 记录每次迭代中最好的适应值 
record_pop_ave = zeros(ger,1);       % 记录种群适应值的平均值





all_one_rank = 0;%初始化为0
while iter <= ger 
    
    disp(['迭代次数:' ,num2str(iter)]);
    
    
    %用一个临时种群去操作 保存struct_pop_new  用 struct_pops去操作
    struct_pops_new = struct_pops;
    
    %第一个迭代时  随机选择一个 
    if iter == 1
        struct_best_indivi = struct_pops;
    else
%         disp(leader(1,(1:1:sizepop)));
        for k=1:sizepop
            struct_best_indivi(k) = struct_pops_archive(leader(1,k));
        end 
    end
    
    %进行特殊处理
    for j=1:sizepop
        %进行特殊的处理
        num_swap = 2;
        rand_index_swap_best = randperm(N,num_swap);
        %从最好的个体那里获取一部分节点
        for z=1:num_swap
            struct_pops(j).per(:,rand_index_swap_best(1,z)) = struct_best_indivi(j).per(:,rand_index_swap_best(1,z));
        end
    end

    %更新个体
    for i=1:sizepop 
         if rand(1,1) > p               %局部搜索
            b = 1-(1-(((ger-iter)/ger).^2)).^0.5;
            index_rand = randperm(sizepop,2);%随机选择两个
            epsilon = b*rand(2,N);%相当于是变化的步长吧
            %%得到新的中间体
            struct_temp_indivi_lo = struct_pops(i);
            index_N_serson = randperm(N,N);%所有结点进行处理
            for j=1:N
                struct_temp_indivi_lo.per(:,index_N_serson(1,j)) = struct_pops(i).per(:,index_N_serson(1,j)) + epsilon(:,j) .*(struct_pops(index_rand(1,1)).per(:,index_N_serson(1,j))- struct_pops(index_rand(1,2)).per(:,index_N_serson(1,j)));

                %同时进行越界处理  2表示x,y坐标
                for k=1:2
                     if struct_temp_indivi_lo.per(k,index_N_serson(1,j)) < pos_limit(1,1) || struct_temp_indivi_lo.per(k,index_N_serson(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_lo.per(k,index_N_serson(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %处理障碍物
                %左上
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(5,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(6,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(8,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4);
                     end
                 end

                %因为上下处理后有一部分会落到中间那部分  所以最后处理中间那个部分
                 %处理矩形上 
                 %障碍物处理  和周围处理  且根据比例来进行分配  上 下各为4  左右各为1 共10
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=35) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=35&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=40)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %上
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 40+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %下
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 35-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==8  %左
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 15-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==9%右
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 35+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end

                 %处理矩形下
                 %障碍物处理  和周围处理  且根据比例来进行分配  上 下各为4  左右各为1 共10
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=35) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=10&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=15)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %上
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 15+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %下
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 10-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==8  %左
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 15-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==9%右
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 35+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end

                 %处理矩形中
                 %障碍物处理  和周围处理  且根据比例来进行分配  只需要分配到左右两侧就行
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=22&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=28) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=35)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3||rem(j,10)==4 %左
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 22-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==5||rem(j,10)==6||rem(j,10)==7||rem(j,10)==8||rem(j,10)==9 %右
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 28+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end
            end
            
            lo_swap_num = 10;
            lo_serson_index = randperm(N,lo_swap_num);
            for k=1:lo_swap_num
                struct_temp_indivi_lo.per(:,lo_serson_index(1,k)) = struct_best_indivi(i).per(:,lo_serson_index(1,k));
            end
            %根据是否支配来判断是否比上一代好 采用精英策略故不用比较
            struct_pops(i) = struct_temp_indivi_lo;%进行替换
         else%进行全局的搜索
            L1 = Levy2(dimension,iter,ger);%得到步长
            index_rand = randperm(N,N);%随机选择两个
            %%得到新的中间体
            struct_temp_indivi_glo = struct_pops(i);
            for j=1:N
                struct_temp_indivi_glo.per(:,index_rand(1,j)) = struct_pops(i).per(:,index_rand(1,j)) + L1' .*(struct_pops(i).per(:,index_rand(1,j))- struct_best_indivi(i).per(:,index_rand(1,j)));%得到临时花朵
            
                %同时进行越界处理  2表示x,y坐标
                for k=1:2
                     if struct_temp_indivi_glo.per(k,index_rand(1,j)) < pos_limit(1,1) || struct_temp_indivi_glo.per(k,index_rand(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_glo.per(k,index_rand(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                
                %处理障碍物
                %左上
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1) <= struct_temp_indivi_glo.per(2,index_rand(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(5,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(6,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(8,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4);
                     end
                 end




                %因为上下处理后有一部分会落到中间那部分  所以最后处理中间那个部分
                 %处理矩形上 
                 %障碍物处理  和周围处理  且根据比例来进行分配  上 下各为4  左右各为1 共10
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=35) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=35&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=40)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %上
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 40+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %下
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 35-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==8  %左
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 15-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==9%右
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 35+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end

                 %处理矩形下
                 %障碍物处理  和周围处理  且根据比例来进行分配  上 下各为4  左右各为1 共10
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=35) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=10&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=15)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %上
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 15+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %下
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 10-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==8  %左
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 15-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==9%右
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 35+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end

                 %处理矩形中
                 %障碍物处理  和周围处理  且根据比例来进行分配  只需要分配到左右两侧就行
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=22&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=28) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=35)%说明在障碍物里面了
                %也是要均匀的放置 均匀放在障碍物的周围 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3||rem(j,10)==4 %左
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 22-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==5||rem(j,10)==6||rem(j,10)==7||rem(j,10)==8||rem(j,10)==9 %右
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 28+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end
            end
            
            %用最好的个体进行下处理
            glo_swap_num = 10;
            glo_serson_index = randperm(N,glo_swap_num);
            for k=1:glo_swap_num
                struct_temp_indivi_glo.per(:,glo_serson_index(1,k)) = struct_best_indivi(i).per(:,glo_serson_index(1,k));
            end
            
           

            %根据是否支配来判断是否比上一代好 采用精英策略故不用比较
            struct_pops(i) = struct_temp_indivi_glo;%进行替换
         end
    end
       
        

	    %整理成大种群
        %进行选择操作
        %整理成大种群
        struct_pops_new_origi((1:1:sizepop)) = struct_pops;%原有的
        struct_pops_new_origi((sizepop+1:1:2*sizepop)) = struct_pops_new;%现在的
        
        
        
        %求各种适应值   记得我们这里求最大的  所以得化成求最大的
        for k=1:2*sizepop
            sensor_mat(1,:) = struct_pops_new_origi(k).per(1,:);
            sensor_mat(2,:) = struct_pops_new_origi(k).per(2,:);
            [origi_new_ifa_cover_fitness(k,1), origi_new_ifa_waste_fitness(k,1)] = get_Grid_cover_unit_and_rate_waste(sensor_mat,struct_pops_new_origi(k).radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,origi_new_ifa_energy_fitness(k,1)] = get_energy_consume(struct_first_init.per,struct_pops_new_origi(k).per,struct_pops_new_origi(k).radius,struct_pops_new_origi(k).energy_init);
        end
        
        
        % 第一行存个体的序号  第二行存被支配的次数，第三行为是否标记，第四行为等级 第5 第6存f1覆盖率 f2节点浪费率适应值 第7行存f3能耗率
        pareto_pop_new_origi = zeros(1 + 6,2 * sizepop);%用来存2倍的种群 第一行存个体序号
        pareto_pop_new_origi(1,(1:1:2*sizepop)) = (1:1:2*sizepop);%节点序号的处理
        

        
        %进行非支配操作
        
        struct_layer_pop = struct('per',[],'nums',[],'layer_num',[]);%结构体类型,个体，同层个数，第二层数
        struct_layer_pops =  repmat(struct_layer_pop,[2*sizepop 1]);%结构体类型的一个种群

        %初始化结构体种群
        %转成行 前1行为节点序号， 第2行为拥挤距离，第3行为等级 第4覆盖率 5浪费率  6能量损耗率
        per_init = zeros(1 + 2 + 3,2*sizepop);
        per_init((1),:) = 0;%初始为0
        per_init(1 + 1,:) = 0;%拥挤距离初始为0
        per_init(1 + 2,:) = 0;%等级初始为0
        per_init((4:1:6),:) = 0;%适应值初始为0
        nums_init = 0;%同一层的个数
        layer_num_init  = 0;%层的序号初始化为0
       

        %组合成结构体
        for i=1:2*sizepop
            struct_layer_pops(i).per  = per_init;
            struct_layer_pops(i).nums = nums_init ;
            struct_layer_pops(i).layer_num = layer_num_init ;

            %顺便把适应值给复制进去
            pareto_pop_new_origi(5,i) = origi_new_ifa_cover_fitness(i,1);
            pareto_pop_new_origi(6,i) = origi_new_ifa_waste_fitness(i,1);
            pareto_pop_new_origi(7,i) = origi_new_ifa_energy_fitness(i,1);
        end
        
         %完成求适应值 并求得对应的最大值和最小值
        f1_fitness_max = max(pareto_pop_new_origi(5,:));
        f1_fitness_min = min(pareto_pop_new_origi(5,:));
        f2_fitness_max = max(pareto_pop_new_origi(6,:));
        f2_fitness_min = min(pareto_pop_new_origi(6,:));
        f3_fitness_max = max(pareto_pop_new_origi(7,:));
        f3_fitness_min = min(pareto_pop_new_origi(7,:));
        
        
        %进行标记等级
        rank = 0 ;%等级标记
        while ~all(pareto_pop_new_origi(3,:)) == 1 
            rank = rank + 1;%下一个等级
            for i=1:2*sizepop
                if pareto_pop_new_origi(3,i) == 1%已经标记了  就不需要进行比较了
                    continue;
                end
                for j=1:2*sizepop
                    if pareto_pop_new_origi(3,j) == 1%已经标记了  就不需要进行比较了
                        continue;
                    end
                    if i==j
                        continue;%如果是自己  那就不比较了
                    end
                    %统计被支配了多少次
                    if pareto_pop_new_origi(5,i) == pareto_pop_new_origi(5,j) && pareto_pop_new_origi(6,i) == pareto_pop_new_origi(6,j) && pareto_pop_new_origi(7,i) == pareto_pop_new_origi(7,j)
                        continue;%如果相等的话  那就没有支配与不支配关系
                    end
                    %统计被支配的次数
                    if pareto_pop_new_origi(5,i) <= pareto_pop_new_origi(5,j) && pareto_pop_new_origi(6,i) >= pareto_pop_new_origi(6,j) && pareto_pop_new_origi(7,i) >= pareto_pop_new_origi(7,j) 
                        pareto_pop_new_origi(2,i) =  pareto_pop_new_origi(2,i) + 1;
                    end
                end    
            end


            %接下来进行  标记和  安排等级
            [~,col] = find(pareto_pop_new_origi(2,:) == 0);%得到row才是有用的
            for k=1:length(col)
                pareto_pop_new_origi(3,col(1,k)) = 1;%标记
                pareto_pop_new_origi(4,col(1,k)) = rank;%等级
                pareto_pop_new_origi(2,col(1,k)) = inf;%仅仅用于标记
            end


            %被支配次数统计进行清零  和 处理 进行下一次统计
            for k=1:2*sizepop
                if pareto_pop_new_origi(2,k) ~= inf
                    pareto_pop_new_origi(2,k) = 0;
                end
            end
        end

        %根据第4列进行排序  行也会跟着动
        pareto_pop_tran = pareto_pop_new_origi';%进行转置
        pareto_pop_sort_new_origin = (sortrows(pareto_pop_tran,4))';%4表示第四列升序排序
        
        %如果标记的等级全是为1，即互相不支配，那么此时是pareto最优 实际用处并不大 还是使用迭代作为循环结束标志
        if all(pareto_pop_sort_new_origin(4,(1:1:sizepop))==1) == 1
            all_one_rank = 1;%表示等级都为1
        end
        
        disp('等级');
        disp(pareto_pop_sort_new_origin(4,(1:1:sizepop)));
        %pause(5);

        %进行分层
        layer_index = 1;%层数从1开始
        count = 1;%个数从0开始
        for i=1:2*sizepop
            %如果是第一次 特殊处理
            if i == 1
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %加入三个函数的适应值
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
                struct_layer_pops(layer_index).layer_num = layer_index;
                continue;
            end
            %不是第一个个体时
            if pareto_pop_sort_new_origin(4,i) == pareto_pop_sort_new_origin(4,i-1)%和上一个个体是同一层
                count = count + 1;
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %加入适应值
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
            else%如果不同层
                %先把上一层的层数、个数统计一下
                struct_layer_pops(layer_index).layer_num = layer_index;
                struct_layer_pops(layer_index).nums = count;

                %重新进行统计
                count = 1;%重新统计
                layer_index = layer_index +  1;%加一层

                %对应的层次和保存个体
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %加入适应值
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
                struct_layer_pops(layer_index).layer_num = layer_index;
            end 
        end
        %最后一层的个数  和 层数少了  补上
         struct_layer_pops(layer_index).nums = count;
         struct_layer_pops(layer_index).layer_num = layer_index;


         %对同一层个体进行处理 补上距离和等级
         layer_max = layer_index;%最大的层数
         for i=1:layer_max
             %根据第一个适应值进行从小到大排序
             temp_mat_f1_tran = struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums))';%转成转置矩阵
             temp_mat_f1_tran_sort = sortrows(temp_mat_f1_tran,4);%根据拥挤距离进行逆序排序
             struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums)) = temp_mat_f1_tran_sort';
             
             
             for j=1:struct_layer_pops(i).nums  %配合下面的公式 
                if j == 1 || j == struct_layer_pops(i).nums
                    struct_layer_pops(i).per(3,j) = i;
                    struct_layer_pops(i).per(2,j) = inf;
                    continue;
                else
                    
                    max_f1_1 = max(max(struct_layer_pops(i).per(4,j-1),struct_layer_pops(i).per(4,j)), struct_layer_pops(i).per(4,j+1));
                    min_f1_1 = min(min(struct_layer_pops(i).per(4,j-1),struct_layer_pops(i).per(4,j)), struct_layer_pops(i).per(4,j+1));
                    max_f2_1 = max(max(struct_layer_pops(i).per(5,j-1),struct_layer_pops(i).per(5,j)), struct_layer_pops(i).per(5,j+1));
                    min_f2_1 = min(min(struct_layer_pops(i).per(5,j-1),struct_layer_pops(i).per(5,j)), struct_layer_pops(i).per(5,j+1));
                    max_f3_1 = max(max(struct_layer_pops(i).per(6,j-1),struct_layer_pops(i).per(6,j)), struct_layer_pops(i).per(6,j+1));
                    min_f3_1 = min(min(struct_layer_pops(i).per(6,j-1),struct_layer_pops(i).per(6,j)), struct_layer_pops(i).per(6,j+1));    
                    
                    struct_layer_pops(i).per(2,j) = (max_f1_1 - min_f1_1) * (max_f2_1 - min_f2_1) * (max_f3_1 - min_f3_1);
                    struct_layer_pops(i).per(3,j) = i;
                end
               

             end
         end
         
         
         
     
     
         
         %进行结构体的处理 对于同一层当中，只需要对距离进行排序 且是从大到小排序
         for i=1:layer_max
             if struct_layer_pops(i).nums == 1
                 continue;%一层只有一个数的话那就没啥比较的了
             else%虽然第2、3、4、5、6行并未使用 但是便于计算，仍使用了它
                 temp_mat_tran = struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums))';%转成转置矩阵
                 temp_mat_sort = sortrows(temp_mat_tran,-2);%根据拥挤距离进行逆序排序
                 temp_mat = temp_mat_sort';%再转置回去
                 struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums)) = temp_mat;%再放回去
             end
         end
         
         
        
         new_pop_index_fitness = zeros(4,sizepop);%存储 第一行序号 2-4对应的适应值


         %存放到新的矩阵中
         indivi_index = 0;%代表的是下标（个体的序号）
         for i=1:layer_max
             for j=1:struct_layer_pops(i).nums
                indivi_index = indivi_index + 1;
                %存到新的矩阵中
                new_pop_index_fitness((1),indivi_index) = struct_layer_pops(i).per((1),j);%序号
                new_pop_index_fitness((2),indivi_index) = struct_layer_pops(i).per((4),j);%个体
                new_pop_index_fitness((3),indivi_index) = struct_layer_pops(i).per((5),j);%个体
                new_pop_index_fitness((4),indivi_index) = struct_layer_pops(i).per((6),j);%个体
                %只存前sizepop个个体
                if indivi_index == sizepop
                    break;
                end
             end
             if indivi_index == sizepop
                break;
             end
         end		
         
     %适应值只为了画图
     figure(12)
     plot3(new_pop_index_fitness(2,:),new_pop_index_fitness(3,:),new_pop_index_fitness(4,:),'.','MarkerSize',20,'color','b');
     grid on;
     pause(0.1);

     %保存最优的sizepop个 到下一次迭代的种群
     for i=1:sizepop
         struct_pops(i) = struct_pops_new_origi(new_pop_index_fitness(1,i));
     end

    %计算每代的平均适应值
    aver_fitness_ger(1,iter) = sum(new_pop_index_fitness(2,:))/sizepop;
    aver_fitness_ger(2,iter) = sum(new_pop_index_fitness(3,:))/sizepop;
    aver_fitness_ger(3,iter) = sum(new_pop_index_fitness(4,:))/sizepop;

    
    
    %解决外部解集的问题
    %把排名靠前的新的sizepop存储到wolf_pop中
    
    %等级为1的个体都是好的个体 加入外部archive
     archive_index = archive_num;%保存上一次的数目
     for j=1:struct_layer_pops(1).nums  
         archive(1,j+archive_index) = struct_layer_pops(1).per(1,j);
         archive(5,j+archive_index) = struct_layer_pops(1).per(4,j);
         archive(6,j+archive_index) = struct_layer_pops(1).per(5,j);
         archive(7,j+archive_index) = struct_layer_pops(1).per(6,j);
         %把个体也放进去 
         struct_pops_archive(:,j+archive_index) = struct_pops_new_origi(:,struct_layer_pops(1).per(1,j));
     end
     archive_num = j+archive_index;%数目
     %序号以struct_archive里面序号为准
     archive(1,(1:1:archive_num)) = (1:1:archive_num);
     
    
     
     
    %进行非支配排序
    for i=1:archive_num 
        for j=1:archive_num 
            if i==j
                continue;%如果是自己  那就不比较了
            end
            %统计被支配了多少次
            if archive(5,i) == archive(5,j) && archive(6,i) == archive(6,j) && archive(7,i) == archive(7,j)
                continue;%如果相等的话  那就没有支配与不支配关系
            end
            %统计被支配的次数
            if archive(5,i) <= archive(5,j) && archive(6,i) >= archive(6,j) && archive(7,i) >= archive(7,j)
                archive(2,i) =  archive(2,i) + 1;
            end
        end    
    end


    %接下来进行  标记和  安排等级
    [row,col] = find(archive(2,(1:1:archive_num)) == 0);%得到row才是有用的
    for k=1:length(col)
        archive(2,col(1,k)) = inf;%仅仅用于标记
    end

%     disp(archive);
    %得到里面的最新的非支配解集    
    archive_temp =  zeros(1+6,2*sizepop);%和pareto设置是一样的
    for k=1:length(col)
        archive_temp(:,k) = archive(:,col(1,k));
        %保存个体
        struct_pops_archive_temp(:,k) = struct_pops_archive(:,col(1,k));
    end
    
    archive_temp_num  = k;%得到当前archive外部的个数
    %重新编序号
    archive_temp(1,(1:1:archive_temp_num)) = (1:1:archive_temp_num);  %以struct_pops_archive_temp的序号为准
    
    
    %对archive_temp进行拥挤距离计算
    %根据其中一个适应值排序
    temp_mat_tran_archive = archive_temp';%转成转置矩阵
    temp_mat_tran_archive_sort = sortrows(temp_mat_tran_archive(1:1:archive_temp_num,:),5);%根据其中一个适应值排序

    archive_temp = temp_mat_tran_archive_sort';

    %求拥挤距离  第四行用于保存拥挤距离
    %第一个和最后一个肯定是会保存下去的 不能让他们太大  合适的数字 让他们处于平均水平吧
    crow_dis_sum = 0;
    for i=1:archive_temp_num
        if i == 1 || i == archive_temp_num
            archive_temp(4,i) = inf;%用数字来替代无穷大  便于计算 往后看
            continue;
        else     
            max_f1 = max(max(archive(5,i),archive(5,i-1)),archive(5,i+1));
            min_f1 = min(min(archive(5,i),archive(5,i-1)),archive(5,i+1));
            max_f2 = max(max(archive(6,i),archive(6,i-1)),archive(6,i+1));
            min_f2 = min(min(archive(6,i),archive(6,i-1)),archive(6,i+1));
            max_f3 = max(max(archive(7,i),archive(7,i-1)),archive(7,i+1));
            min_f3 = min(min(archive(7,i),archive(7,i-1)),archive(7,i+1));
            %拥挤距离
            archive_temp(4,i) = (max_f1 - min_f1 + 0.0001) * (max_f2 - min_f2 + 0.0001) * (max_f3 - min_f3 + 0.0001);
            crow_dis_sum = crow_dis_sum + archive_temp(4,i);
        end
    end
    %解决第一个 和最后一个  为了计算
    archive_temp(4,1) = crow_dis_sum/(archive_temp_num-2);
    archive_temp(4,archive_temp_num) = archive_temp(4,1);
    
    
    %再根据拥挤距离进行逆序排序
    temp_mat_tran_dis_archive = archive_temp';%转成转置矩阵
    temp_mat_tran_archive_dis_sort = sortrows(temp_mat_tran_dis_archive,-4);%根据其中一个适应值排序
    archive_temp = temp_mat_tran_archive_dis_sort';
   
    

    
    archive((1:1:7),(1:1:2*sizepop)) =  0;%清0
    %如果超过了限制的个数   则得把最优的前archive_num_limit保留下来
    if archive_temp_num > archive_num_limit
        %把排好序的前archive_num_limit放入到archive中
        for i=1:archive_num_limit
            archive(:,i) = archive_temp(:,i);
            %保存个体
            struct_pops_archive(:,i) = struct_pops_archive_temp(:,i);
        end
        %标记当前archive中个体的数目
        archive_num = archive_num_limit;
        %序号重新编
        archive(1,(1:1:archive_num)) = (1:1:archive_num);
    else
        for i=1:archive_temp_num
            archive(:,i) = archive_temp(:,i);
            %保存个体
            struct_pops_archive(:,i) = struct_pops_archive_temp(:,i);
        end
        archive_num = archive_temp_num;
        %序号重新编
        archive(1,(1:1:archive_num)) = (1:1:archive_num);
    end
    disp(['archive个数:' ,num2str(archive_num)]);
    
    figure(11);
    plot3(archive(5,(1:1:archive_num)),archive(6,(1:1:archive_num)),archive(7,(1:1:archive_num)),'.','MarkerSize',20,'color','b');
    grid on;
    pause(0.01);
    
    %把得到的archive 根据拥挤距离轮盘赌法  初始化leader
    %先构造一个轮盘赌法
    totalfit = sum(archive(4,(1:1:archive_num)));%得到总的适应度值
    p_fitvalue = archive(4,(1:1:archive_num))'./totalfit;%得到每个适应度值占总的适应度的比例

    
    %得到概率带
    p_sum = 0;%初始为0
    p_fitvalue_ban = zeros(archive_num+1,1);%概率带预分配内存  还得考虑第一元素为0
    p_fitvalue_ban(1,1) = 0;%第一个元素
    for i=1:archive_num
        p_sum = p_sum + p_fitvalue(i,1);
        p_fitvalue_ban((i+1),1) = p_sum;
    end

    space_count = zeros(archive_num,1);%用于统计每一个间隔内的个数
    
    for i=1:archive_num_limit%sizepop这里同时代表的是概率带的间隔数目
        rand_p = rand(1,1);
        for j=1:archive_num%sizepop这里同时代表的是概率带的间隔数目
            if rand_p >= p_fitvalue_ban(j,1) && rand_p < p_fitvalue_ban(j+1,1)
                space_count(j,1) = space_count(j,1)+1;
                break;%停止继续搜索
            end
        end
    end

%     这就是选择
%   创建一个新的种群并且把相应的适应值得到
    leader = zeros(1+6,sizepop);%和pareto设置是一样的  4行是拥挤距离
    k = 1;
    for j=1:archive_num
        while space_count(j,1)>0
            leader(:,k) = archive(:,j);%得到相应的新种群
            k = k+1;%下标走动
            space_count(j,1) = space_count(j,1)-1;%数目减少
        end
    end
    
    disp(leader(1,:));
    %把archive里面的除了数据(序号、三个适应值)之外的全部清0 1 5 6 7不需要
    archive((2:1:4),(1:1:2*sizepop)) =  0;%清0
    
    iter = iter+1; 
end

%画三个适应值的图
figure(110);
plot((1:1:ger),aver_fitness_ger(1,:));
hold on;

figure(111);
plot((1:1:ger),aver_fitness_ger(2,:));
hold on;

figure(112);
plot((1:1:ger),aver_fitness_ger(3,:));
hold on;

%保存两重要数据
save ifa_aver_fitness_ger.mat aver_fitness_ger;%三个平均值
save ifa_WSN_fitness.mat new_pop_index_fitness;%最后的结果
