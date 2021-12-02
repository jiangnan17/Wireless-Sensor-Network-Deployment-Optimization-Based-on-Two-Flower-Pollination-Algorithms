%%������
clc;
clear ;
close all;
%ɾ����Ӧ���ļ�

global N;
global M;
global L;
global W;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global ger;

p=0.8;%�ж��Ƿ���ȫ���Ż����Ǿֲ��Ż�

L = 50;%��
W = 50;%��
%����1ƽ����һ������
M = 2500;%��������
r_max = 7;%��֪�뾶Ϊ5
r_mid = 6;
r_min = 5;
energy_max = 100;%��������
energy_mid = 90;
energy_min = 80;

per_sersons_radius_type = [r_max,r_mid,r_min];
%�������Ϊ5��ʣ��ΪС
N = 25;%30���������ڵ�
sizepop = 200;%��Ⱥ��ģ
dimension = 2;% �ռ�ά��  ǰ�з�x��y�������зŰ뾶
ger = 2;% ����������
pos_limit = [0, 50];            % ����λ�ò�������
%��������
r_max_num = 1;%���Ϊ1-5
r_mid_num = 2;%���Ϊ6-10
r_min_num = N - r_max_num - r_mid_num;%���Ϊ11-N

aver_fitness_ger = zeros(3,ger);%��������Ӧֵÿ����ƽ��ֵ
struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);%�ṹ������

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%��ʱ��һ����Ⱥ

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);


%������ε��ĸ���
syms x y;%�ȶ���һ������
%���Ͻ�
k1 = 1;
b1 = 35;
x1_up = solve(k1*x+b1==50,x);%���Ͻǵ�б�ߵ��ϸ�����
y1_down = solve(k1*0+b1==y,y);

%���½�
k2 = -1;
b2 = 15;
y2_up = solve(k2*0+b2==y,y);
x2_down = solve(k2*x+b2==0,x);


%���Ͻ�
k3 = -1;
b3 = 85;
x3_up = solve(k3*x+b3==50,x);
y3_down = solve(k3*50+b3==y,y);

%���½�
k4 = 1;
b4 = -35;
y4_up = solve(k4*50+b4==y,y);
x4_down = solve(k4*x+b4==0,x);

%����������֤��ϣ���ȫ��ȷ
point = zeros(8,2);%�洢��Щ��  ����  ��������
point(1,:) = [x1_up,50];
point(2,:) = [0,y1_down];
point(3,:) = [0,y2_up];
point(4,:) = [x2_down,0];
point(5,:) = [x3_up,50];
point(6,:) = [50,y3_down];
point(7,:) = [50,y4_up];
point(8,:) = [x4_down,0];

load struct_pop_public.mat;%���ظ���Ⱥ
struct_pops = struct_pop_public;%�õ���Ⱥ����
struct_pops_new_origi = repmat(struct_pop_per,[1 2*sizepop]);%����ԭ����Ⱥ���µ���Ⱥ

struct_pops_archive = struct_pops_new_origi;%���ڱ���archive����ĸ���


struct_pops_archive_temp = struct_pops_new_origi;%���ڱ���archive_temp����ĸ���


load struct_first_init_public.mat%�����ʼ��һ����������
struct_first_init = struct_first_init_public;%�õ���ʼ����������


%%��ʼ�Ĳ����ͼ  �õ�һ��������ȥ��ʼ��ͼ

%��������������
X_mat = (0:1:50);%x����
Y_mat = (0:1:50);%y����
Grid_cen_x = zeros(1,L/1);%�������ĵ�x����
Grid_cen_y = zeros(1,W/1);%�������ĵ�y����
%ǰ���������֮�ͳ���2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end


%%�Ѻ������궪��һ����ά������
%����ת����  ��һ�з�x�� �ڶ��з�y�ᣬͬһ�з�һ��������
%���ȴ濿��x����ĵ�һ�У�Ȼ�����ϴ�ڶ���
%������������
Grid_cen_x_and_y = zeros(L,W,2);%��2500���������ģ�����ÿ������������x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1����x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%��y����ŵ��ڶ���
    end
end


x_pos = struct_first_init.per(1,:);%��һ������   ���Ӽ��ǽ�  ��x����
y_pos = struct_first_init.per(2,:);%��һ������   ���Ӽ��ǽ�  ��y����
sersors_r = struct_first_init.radius;%��һ������ 




%���л�ͼ
% figure(1);
% draw_circle(x_pos,y_pos,sersors_r);
% title('��ʼ������ͼ');
% hold on;





%���������ĵ������ŵ�������
sensor_mat = zeros(2,N);%Ԥ�����ڴ�
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

%������ֽڵ����Ŀ
per_sersons_num = struct_pops(1).sersons_num;
%������ֽڵ�İ뾶


%%��ʼ���õ����ϸ��ʺͽڵ��˷���
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['��ʼ���ĸ����ʣ�',num2str(cover_rate)]);
disp(['��ʼ�����˷��ʣ�',num2str(waste_rate)]);



%%������ͨ��  ��һ���ڵ��ʼ��ʱ����ͨ��
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('��ͨ');
else
    disp('����ͨ');
end




%%��ʼ����Ⱥ��ʷֵΪ  ����С   infΪ�����

best_fitness = -inf;                         % ��Ⱥ��ʷ�����Ӧ��  
struct_best_indivi = struct_pop_per;                 % �����������
struct_best_indivi_fitness = struct('cover_rate',[],'waste_rate',[],'energy_rate',[],'function_rate',[]);%��Ӧֵ�ṹ������
struct_best_indivi_fitness_all = repmat(struct_best_indivi_fitness,[1 1]);%Ԥ�����ڴ�
%%�������Ѿ���ͼ��
%%����Ϊ������ʼ��ʱ��50�����ӿ�ʼ��λ��


archive = zeros(1+6,2*sizepop);%��pareto������һ����  4����ӵ������
archive_num = 0;%����ʵ�ڵķ�֧��
archive_num_limit = sizepop;%����쵼�ߵ�����
%����������޷�����  ����һ���ܴ������ inf_value  ������̶ķ����������


origi_new_ifa_cover_fitness = zeros(2*sizepop,1);%���帲����
origi_new_ifa_waste_fitness = zeros(2*sizepop,1);%�����˷���
origi_new_ifa_energy_fitness = zeros(2*sizepop,1);%�����ܺ���  ����1/xx��ʽ 

ifa_function_fitness = zeros(2*sizepop,1);%���������Ľ��
ifa_fitness_temp = zeros(2*sizepop,1);%��ʱ���浱ǰ��Ӧ��
%% Ⱥ�����
iter = 1;
record_ger = zeros(ger, 1);          % ��¼ÿ�ε�������õ���Ӧֵ 
record_pop_ave = zeros(ger,1);       % ��¼��Ⱥ��Ӧֵ��ƽ��ֵ





all_one_rank = 0;%��ʼ��Ϊ0
while iter <= ger 
    
    disp(['��������:' ,num2str(iter)]);
    
    
    %��һ����ʱ��Ⱥȥ���� ����struct_pop_new  �� struct_popsȥ����
    struct_pops_new = struct_pops;
    
    %��һ������ʱ  ���ѡ��һ�� 
    if iter == 1
        struct_best_indivi = struct_pops;
    else
%         disp(leader(1,(1:1:sizepop)));
        for k=1:sizepop
            struct_best_indivi(k) = struct_pops_archive(leader(1,k));
        end 
    end
    
    %�������⴦��
    for j=1:sizepop
        %��������Ĵ���
        num_swap = 2;
        rand_index_swap_best = randperm(N,num_swap);
        %����õĸ��������ȡһ���ֽڵ�
        for z=1:num_swap
            struct_pops(j).per(:,rand_index_swap_best(1,z)) = struct_best_indivi(j).per(:,rand_index_swap_best(1,z));
        end
    end

    %���¸���
    for i=1:sizepop 
         if rand(1,1) > p               %�ֲ�����
            b = 1-(1-(((ger-iter)/ger).^2)).^0.5;
            index_rand = randperm(sizepop,2);%���ѡ������
            epsilon = b*rand(2,N);%�൱���Ǳ仯�Ĳ�����
            %%�õ��µ��м���
            struct_temp_indivi_lo = struct_pops(i);
            index_N_serson = randperm(N,N);%���н����д���
            for j=1:N
                struct_temp_indivi_lo.per(:,index_N_serson(1,j)) = struct_pops(i).per(:,index_N_serson(1,j)) + epsilon(:,j) .*(struct_pops(index_rand(1,1)).per(:,index_N_serson(1,j))- struct_pops(index_rand(1,2)).per(:,index_N_serson(1,j)));

                %ͬʱ����Խ�紦��  2��ʾx,y����
                for k=1:2
                     if struct_temp_indivi_lo.per(k,index_N_serson(1,j)) < pos_limit(1,1) || struct_temp_indivi_lo.per(k,index_N_serson(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_lo.per(k,index_N_serson(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %�����ϰ���
                %����
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(5,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(6,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(8,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4);
                     end
                 end

                %��Ϊ���´������һ���ֻ��䵽�м��ǲ���  ����������м��Ǹ�����
                 %��������� 
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  �� �¸�Ϊ4  ���Ҹ�Ϊ1 ��10
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=35) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=35&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=40)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %��
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 40+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %��
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 35-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==8  %��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 15-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==9%��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 35+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end

                 %���������
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  �� �¸�Ϊ4  ���Ҹ�Ϊ1 ��10
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=35) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=10&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=15)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %��
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 15+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %��
                        struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = 10-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==8  %��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 15-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==9%��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 35+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end

                 %���������
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  ֻ��Ҫ���䵽�����������
                if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=22&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=28) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=15&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=35)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3||rem(j,10)==4 %��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 22-rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    elseif rem(j,10)==5||rem(j,10)==6||rem(j,10)==7||rem(j,10)==8||rem(j,10)==9 %��
                        struct_temp_indivi_lo.per(1,index_N_serson(1,j)) = 28+rand(1,1)*struct_temp_indivi_lo.radius(1,index_N_serson(1,j));
                    end
                end
            end
            
            lo_swap_num = 10;
            lo_serson_index = randperm(N,lo_swap_num);
            for k=1:lo_swap_num
                struct_temp_indivi_lo.per(:,lo_serson_index(1,k)) = struct_best_indivi(i).per(:,lo_serson_index(1,k));
            end
            %�����Ƿ�֧�����ж��Ƿ����һ���� ���þ�Ӣ���Թʲ��ñȽ�
            struct_pops(i) = struct_temp_indivi_lo;%�����滻
         else%����ȫ�ֵ�����
            L1 = Levy2(dimension,iter,ger);%�õ�����
            index_rand = randperm(N,N);%���ѡ������
            %%�õ��µ��м���
            struct_temp_indivi_glo = struct_pops(i);
            for j=1:N
                struct_temp_indivi_glo.per(:,index_rand(1,j)) = struct_pops(i).per(:,index_rand(1,j)) + L1' .*(struct_pops(i).per(:,index_rand(1,j))- struct_best_indivi(i).per(:,index_rand(1,j)));%�õ���ʱ����
            
                %ͬʱ����Խ�紦��  2��ʾx,y����
                for k=1:2
                     if struct_temp_indivi_glo.per(k,index_rand(1,j)) < pos_limit(1,1) || struct_temp_indivi_glo.per(k,index_rand(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_glo.per(k,index_rand(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                
                %�����ϰ���
                %����
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1) <= struct_temp_indivi_glo.per(2,index_rand(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(5,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(6,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(8,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4);
                     end
                 end




                %��Ϊ���´������һ���ֻ��䵽�м��ǲ���  ����������м��Ǹ�����
                 %��������� 
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  �� �¸�Ϊ4  ���Ҹ�Ϊ1 ��10
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=35) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=35&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=40)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %��
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 40+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %��
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 35-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==8  %��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 15-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==9%��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 35+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end

                 %���������
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  �� �¸�Ϊ4  ���Ҹ�Ϊ1 ��10
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=35) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=10&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=15)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3 %��
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 15+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==4||rem(j,10)==5||rem(j,10)==6||rem(j,10)==7 %��
                        struct_temp_indivi_glo.per(2,index_rand(1,j)) = 10-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==8  %��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 15-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==9%��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 35+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end

                 %���������
                 %�ϰ��ﴦ��  ����Χ����  �Ҹ��ݱ��������з���  ֻ��Ҫ���䵽�����������
                if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=22&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=28) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=15&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=35)%˵�����ϰ���������
                %Ҳ��Ҫ���ȵķ��� ���ȷ����ϰ������Χ 
                    if rem(j,10)==0||rem(j,10)==1||rem(j,10)==2||rem(j,10)==3||rem(j,10)==4 %��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 22-rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    elseif rem(j,10)==5||rem(j,10)==6||rem(j,10)==7||rem(j,10)==8||rem(j,10)==9 %��
                        struct_temp_indivi_glo.per(1,index_rand(1,j)) = 28+rand(1,1)*struct_temp_indivi_glo.radius(1,index_rand(1,j));
                    end
                end
            end
            
            %����õĸ�������´���
            glo_swap_num = 10;
            glo_serson_index = randperm(N,glo_swap_num);
            for k=1:glo_swap_num
                struct_temp_indivi_glo.per(:,glo_serson_index(1,k)) = struct_best_indivi(i).per(:,glo_serson_index(1,k));
            end
            
           

            %�����Ƿ�֧�����ж��Ƿ����һ���� ���þ�Ӣ���Թʲ��ñȽ�
            struct_pops(i) = struct_temp_indivi_glo;%�����滻
         end
    end
       
        

	    %����ɴ���Ⱥ
        %����ѡ�����
        %����ɴ���Ⱥ
        struct_pops_new_origi((1:1:sizepop)) = struct_pops;%ԭ�е�
        struct_pops_new_origi((sizepop+1:1:2*sizepop)) = struct_pops_new;%���ڵ�
        
        
        
        %�������Ӧֵ   �ǵ���������������  ���Եû���������
        for k=1:2*sizepop
            sensor_mat(1,:) = struct_pops_new_origi(k).per(1,:);
            sensor_mat(2,:) = struct_pops_new_origi(k).per(2,:);
            [origi_new_ifa_cover_fitness(k,1), origi_new_ifa_waste_fitness(k,1)] = get_Grid_cover_unit_and_rate_waste(sensor_mat,struct_pops_new_origi(k).radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
            [~,origi_new_ifa_energy_fitness(k,1)] = get_energy_consume(struct_first_init.per,struct_pops_new_origi(k).per,struct_pops_new_origi(k).radius,struct_pops_new_origi(k).energy_init);
        end
        
        
        % ��һ�д��������  �ڶ��д汻֧��Ĵ�����������Ϊ�Ƿ��ǣ�������Ϊ�ȼ� ��5 ��6��f1������ f2�ڵ��˷�����Ӧֵ ��7�д�f3�ܺ���
        pareto_pop_new_origi = zeros(1 + 6,2 * sizepop);%������2������Ⱥ ��һ�д�������
        pareto_pop_new_origi(1,(1:1:2*sizepop)) = (1:1:2*sizepop);%�ڵ���ŵĴ���
        

        
        %���з�֧�����
        
        struct_layer_pop = struct('per',[],'nums',[],'layer_num',[]);%�ṹ������,���壬ͬ��������ڶ�����
        struct_layer_pops =  repmat(struct_layer_pop,[2*sizepop 1]);%�ṹ�����͵�һ����Ⱥ

        %��ʼ���ṹ����Ⱥ
        %ת���� ǰ1��Ϊ�ڵ���ţ� ��2��Ϊӵ�����룬��3��Ϊ�ȼ� ��4������ 5�˷���  6���������
        per_init = zeros(1 + 2 + 3,2*sizepop);
        per_init((1),:) = 0;%��ʼΪ0
        per_init(1 + 1,:) = 0;%ӵ�������ʼΪ0
        per_init(1 + 2,:) = 0;%�ȼ���ʼΪ0
        per_init((4:1:6),:) = 0;%��Ӧֵ��ʼΪ0
        nums_init = 0;%ͬһ��ĸ���
        layer_num_init  = 0;%�����ų�ʼ��Ϊ0
       

        %��ϳɽṹ��
        for i=1:2*sizepop
            struct_layer_pops(i).per  = per_init;
            struct_layer_pops(i).nums = nums_init ;
            struct_layer_pops(i).layer_num = layer_num_init ;

            %˳�����Ӧֵ�����ƽ�ȥ
            pareto_pop_new_origi(5,i) = origi_new_ifa_cover_fitness(i,1);
            pareto_pop_new_origi(6,i) = origi_new_ifa_waste_fitness(i,1);
            pareto_pop_new_origi(7,i) = origi_new_ifa_energy_fitness(i,1);
        end
        
         %�������Ӧֵ ����ö�Ӧ�����ֵ����Сֵ
        f1_fitness_max = max(pareto_pop_new_origi(5,:));
        f1_fitness_min = min(pareto_pop_new_origi(5,:));
        f2_fitness_max = max(pareto_pop_new_origi(6,:));
        f2_fitness_min = min(pareto_pop_new_origi(6,:));
        f3_fitness_max = max(pareto_pop_new_origi(7,:));
        f3_fitness_min = min(pareto_pop_new_origi(7,:));
        
        
        %���б�ǵȼ�
        rank = 0 ;%�ȼ����
        while ~all(pareto_pop_new_origi(3,:)) == 1 
            rank = rank + 1;%��һ���ȼ�
            for i=1:2*sizepop
                if pareto_pop_new_origi(3,i) == 1%�Ѿ������  �Ͳ���Ҫ���бȽ���
                    continue;
                end
                for j=1:2*sizepop
                    if pareto_pop_new_origi(3,j) == 1%�Ѿ������  �Ͳ���Ҫ���бȽ���
                        continue;
                    end
                    if i==j
                        continue;%������Լ�  �ǾͲ��Ƚ���
                    end
                    %ͳ�Ʊ�֧���˶��ٴ�
                    if pareto_pop_new_origi(5,i) == pareto_pop_new_origi(5,j) && pareto_pop_new_origi(6,i) == pareto_pop_new_origi(6,j) && pareto_pop_new_origi(7,i) == pareto_pop_new_origi(7,j)
                        continue;%�����ȵĻ�  �Ǿ�û��֧���벻֧���ϵ
                    end
                    %ͳ�Ʊ�֧��Ĵ���
                    if pareto_pop_new_origi(5,i) <= pareto_pop_new_origi(5,j) && pareto_pop_new_origi(6,i) >= pareto_pop_new_origi(6,j) && pareto_pop_new_origi(7,i) >= pareto_pop_new_origi(7,j) 
                        pareto_pop_new_origi(2,i) =  pareto_pop_new_origi(2,i) + 1;
                    end
                end    
            end


            %����������  ��Ǻ�  ���ŵȼ�
            [~,col] = find(pareto_pop_new_origi(2,:) == 0);%�õ�row�������õ�
            for k=1:length(col)
                pareto_pop_new_origi(3,col(1,k)) = 1;%���
                pareto_pop_new_origi(4,col(1,k)) = rank;%�ȼ�
                pareto_pop_new_origi(2,col(1,k)) = inf;%�������ڱ��
            end


            %��֧�����ͳ�ƽ�������  �� ���� ������һ��ͳ��
            for k=1:2*sizepop
                if pareto_pop_new_origi(2,k) ~= inf
                    pareto_pop_new_origi(2,k) = 0;
                end
            end
        end

        %���ݵ�4�н�������  ��Ҳ����Ŷ�
        pareto_pop_tran = pareto_pop_new_origi';%����ת��
        pareto_pop_sort_new_origin = (sortrows(pareto_pop_tran,4))';%4��ʾ��������������
        
        %�����ǵĵȼ�ȫ��Ϊ1�������಻֧�䣬��ô��ʱ��pareto���� ʵ���ô������� ����ʹ�õ�����Ϊѭ��������־
        if all(pareto_pop_sort_new_origin(4,(1:1:sizepop))==1) == 1
            all_one_rank = 1;%��ʾ�ȼ���Ϊ1
        end
        
        disp('�ȼ�');
        disp(pareto_pop_sort_new_origin(4,(1:1:sizepop)));
        %pause(5);

        %���зֲ�
        layer_index = 1;%������1��ʼ
        count = 1;%������0��ʼ
        for i=1:2*sizepop
            %����ǵ�һ�� ���⴦��
            if i == 1
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %����������������Ӧֵ
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
                struct_layer_pops(layer_index).layer_num = layer_index;
                continue;
            end
            %���ǵ�һ������ʱ
            if pareto_pop_sort_new_origin(4,i) == pareto_pop_sort_new_origin(4,i-1)%����һ��������ͬһ��
                count = count + 1;
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %������Ӧֵ
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
            else%�����ͬ��
                %�Ȱ���һ��Ĳ���������ͳ��һ��
                struct_layer_pops(layer_index).layer_num = layer_index;
                struct_layer_pops(layer_index).nums = count;

                %���½���ͳ��
                count = 1;%����ͳ��
                layer_index = layer_index +  1;%��һ��

                %��Ӧ�Ĳ�κͱ������
                struct_layer_pops(layer_index).per((1),count) = pareto_pop_sort_new_origin((1),i);
                %������Ӧֵ
                struct_layer_pops(layer_index).per((4),count) = pareto_pop_sort_new_origin((5),i);
                struct_layer_pops(layer_index).per((5),count) = pareto_pop_sort_new_origin((6),i);
                struct_layer_pops(layer_index).per((6),count) = pareto_pop_sort_new_origin((7),i);
                struct_layer_pops(layer_index).layer_num = layer_index;
            end 
        end
        %���һ��ĸ���  �� ��������  ����
         struct_layer_pops(layer_index).nums = count;
         struct_layer_pops(layer_index).layer_num = layer_index;


         %��ͬһ�������д��� ���Ͼ���͵ȼ�
         layer_max = layer_index;%���Ĳ���
         for i=1:layer_max
             %���ݵ�һ����Ӧֵ���д�С��������
             temp_mat_f1_tran = struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums))';%ת��ת�þ���
             temp_mat_f1_tran_sort = sortrows(temp_mat_f1_tran,4);%����ӵ�����������������
             struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums)) = temp_mat_f1_tran_sort';
             
             
             for j=1:struct_layer_pops(i).nums  %�������Ĺ�ʽ 
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
         
         
         
     
     
         
         %���нṹ��Ĵ��� ����ͬһ�㵱�У�ֻ��Ҫ�Ծ���������� ���ǴӴ�С����
         for i=1:layer_max
             if struct_layer_pops(i).nums == 1
                 continue;%һ��ֻ��һ�����Ļ��Ǿ�ûɶ�Ƚϵ���
             else%��Ȼ��2��3��4��5��6�в�δʹ�� ���Ǳ��ڼ��㣬��ʹ������
                 temp_mat_tran = struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums))';%ת��ת�þ���
                 temp_mat_sort = sortrows(temp_mat_tran,-2);%����ӵ�����������������
                 temp_mat = temp_mat_sort';%��ת�û�ȥ
                 struct_layer_pops(i).per((1:1:6),(1:1:struct_layer_pops(i).nums)) = temp_mat;%�ٷŻ�ȥ
             end
         end
         
         
        
         new_pop_index_fitness = zeros(4,sizepop);%�洢 ��һ����� 2-4��Ӧ����Ӧֵ


         %��ŵ��µľ�����
         indivi_index = 0;%��������±꣨�������ţ�
         for i=1:layer_max
             for j=1:struct_layer_pops(i).nums
                indivi_index = indivi_index + 1;
                %�浽�µľ�����
                new_pop_index_fitness((1),indivi_index) = struct_layer_pops(i).per((1),j);%���
                new_pop_index_fitness((2),indivi_index) = struct_layer_pops(i).per((4),j);%����
                new_pop_index_fitness((3),indivi_index) = struct_layer_pops(i).per((5),j);%����
                new_pop_index_fitness((4),indivi_index) = struct_layer_pops(i).per((6),j);%����
                %ֻ��ǰsizepop������
                if indivi_index == sizepop
                    break;
                end
             end
             if indivi_index == sizepop
                break;
             end
         end		
         
     %��ӦֵֻΪ�˻�ͼ
     figure(12)
     plot3(new_pop_index_fitness(2,:),new_pop_index_fitness(3,:),new_pop_index_fitness(4,:),'.','MarkerSize',20,'color','b');
     grid on;
     pause(0.1);

     %�������ŵ�sizepop�� ����һ�ε�������Ⱥ
     for i=1:sizepop
         struct_pops(i) = struct_pops_new_origi(new_pop_index_fitness(1,i));
     end

    %����ÿ����ƽ����Ӧֵ
    aver_fitness_ger(1,iter) = sum(new_pop_index_fitness(2,:))/sizepop;
    aver_fitness_ger(2,iter) = sum(new_pop_index_fitness(3,:))/sizepop;
    aver_fitness_ger(3,iter) = sum(new_pop_index_fitness(4,:))/sizepop;

    
    
    %����ⲿ�⼯������
    %��������ǰ���µ�sizepop�洢��wolf_pop��
    
    %�ȼ�Ϊ1�ĸ��嶼�Ǻõĸ��� �����ⲿarchive
     archive_index = archive_num;%������һ�ε���Ŀ
     for j=1:struct_layer_pops(1).nums  
         archive(1,j+archive_index) = struct_layer_pops(1).per(1,j);
         archive(5,j+archive_index) = struct_layer_pops(1).per(4,j);
         archive(6,j+archive_index) = struct_layer_pops(1).per(5,j);
         archive(7,j+archive_index) = struct_layer_pops(1).per(6,j);
         %�Ѹ���Ҳ�Ž�ȥ 
         struct_pops_archive(:,j+archive_index) = struct_pops_new_origi(:,struct_layer_pops(1).per(1,j));
     end
     archive_num = j+archive_index;%��Ŀ
     %�����struct_archive�������Ϊ׼
     archive(1,(1:1:archive_num)) = (1:1:archive_num);
     
    
     
     
    %���з�֧������
    for i=1:archive_num 
        for j=1:archive_num 
            if i==j
                continue;%������Լ�  �ǾͲ��Ƚ���
            end
            %ͳ�Ʊ�֧���˶��ٴ�
            if archive(5,i) == archive(5,j) && archive(6,i) == archive(6,j) && archive(7,i) == archive(7,j)
                continue;%�����ȵĻ�  �Ǿ�û��֧���벻֧���ϵ
            end
            %ͳ�Ʊ�֧��Ĵ���
            if archive(5,i) <= archive(5,j) && archive(6,i) >= archive(6,j) && archive(7,i) >= archive(7,j)
                archive(2,i) =  archive(2,i) + 1;
            end
        end    
    end


    %����������  ��Ǻ�  ���ŵȼ�
    [row,col] = find(archive(2,(1:1:archive_num)) == 0);%�õ�row�������õ�
    for k=1:length(col)
        archive(2,col(1,k)) = inf;%�������ڱ��
    end

%     disp(archive);
    %�õ���������µķ�֧��⼯    
    archive_temp =  zeros(1+6,2*sizepop);%��pareto������һ����
    for k=1:length(col)
        archive_temp(:,k) = archive(:,col(1,k));
        %�������
        struct_pops_archive_temp(:,k) = struct_pops_archive(:,col(1,k));
    end
    
    archive_temp_num  = k;%�õ���ǰarchive�ⲿ�ĸ���
    %���±����
    archive_temp(1,(1:1:archive_temp_num)) = (1:1:archive_temp_num);  %��struct_pops_archive_temp�����Ϊ׼
    
    
    %��archive_temp����ӵ���������
    %��������һ����Ӧֵ����
    temp_mat_tran_archive = archive_temp';%ת��ת�þ���
    temp_mat_tran_archive_sort = sortrows(temp_mat_tran_archive(1:1:archive_temp_num,:),5);%��������һ����Ӧֵ����

    archive_temp = temp_mat_tran_archive_sort';

    %��ӵ������  ���������ڱ���ӵ������
    %��һ�������һ���϶��ǻᱣ����ȥ�� ����������̫��  ���ʵ����� �����Ǵ���ƽ��ˮƽ��
    crow_dis_sum = 0;
    for i=1:archive_temp_num
        if i == 1 || i == archive_temp_num
            archive_temp(4,i) = inf;%����������������  ���ڼ��� ����
            continue;
        else     
            max_f1 = max(max(archive(5,i),archive(5,i-1)),archive(5,i+1));
            min_f1 = min(min(archive(5,i),archive(5,i-1)),archive(5,i+1));
            max_f2 = max(max(archive(6,i),archive(6,i-1)),archive(6,i+1));
            min_f2 = min(min(archive(6,i),archive(6,i-1)),archive(6,i+1));
            max_f3 = max(max(archive(7,i),archive(7,i-1)),archive(7,i+1));
            min_f3 = min(min(archive(7,i),archive(7,i-1)),archive(7,i+1));
            %ӵ������
            archive_temp(4,i) = (max_f1 - min_f1 + 0.0001) * (max_f2 - min_f2 + 0.0001) * (max_f3 - min_f3 + 0.0001);
            crow_dis_sum = crow_dis_sum + archive_temp(4,i);
        end
    end
    %�����һ�� �����һ��  Ϊ�˼���
    archive_temp(4,1) = crow_dis_sum/(archive_temp_num-2);
    archive_temp(4,archive_temp_num) = archive_temp(4,1);
    
    
    %�ٸ���ӵ�����������������
    temp_mat_tran_dis_archive = archive_temp';%ת��ת�þ���
    temp_mat_tran_archive_dis_sort = sortrows(temp_mat_tran_dis_archive,-4);%��������һ����Ӧֵ����
    archive_temp = temp_mat_tran_archive_dis_sort';
   
    

    
    archive((1:1:7),(1:1:2*sizepop)) =  0;%��0
    %������������Ƶĸ���   ��ð����ŵ�ǰarchive_num_limit��������
    if archive_temp_num > archive_num_limit
        %���ź����ǰarchive_num_limit���뵽archive��
        for i=1:archive_num_limit
            archive(:,i) = archive_temp(:,i);
            %�������
            struct_pops_archive(:,i) = struct_pops_archive_temp(:,i);
        end
        %��ǵ�ǰarchive�и������Ŀ
        archive_num = archive_num_limit;
        %������±�
        archive(1,(1:1:archive_num)) = (1:1:archive_num);
    else
        for i=1:archive_temp_num
            archive(:,i) = archive_temp(:,i);
            %�������
            struct_pops_archive(:,i) = struct_pops_archive_temp(:,i);
        end
        archive_num = archive_temp_num;
        %������±�
        archive(1,(1:1:archive_num)) = (1:1:archive_num);
    end
    disp(['archive����:' ,num2str(archive_num)]);
    
    figure(11);
    plot3(archive(5,(1:1:archive_num)),archive(6,(1:1:archive_num)),archive(7,(1:1:archive_num)),'.','MarkerSize',20,'color','b');
    grid on;
    pause(0.01);
    
    %�ѵõ���archive ����ӵ���������̶ķ�  ��ʼ��leader
    %�ȹ���һ�����̶ķ�
    totalfit = sum(archive(4,(1:1:archive_num)));%�õ��ܵ���Ӧ��ֵ
    p_fitvalue = archive(4,(1:1:archive_num))'./totalfit;%�õ�ÿ����Ӧ��ֵռ�ܵ���Ӧ�ȵı���

    
    %�õ����ʴ�
    p_sum = 0;%��ʼΪ0
    p_fitvalue_ban = zeros(archive_num+1,1);%���ʴ�Ԥ�����ڴ�  ���ÿ��ǵ�һԪ��Ϊ0
    p_fitvalue_ban(1,1) = 0;%��һ��Ԫ��
    for i=1:archive_num
        p_sum = p_sum + p_fitvalue(i,1);
        p_fitvalue_ban((i+1),1) = p_sum;
    end

    space_count = zeros(archive_num,1);%����ͳ��ÿһ������ڵĸ���
    
    for i=1:archive_num_limit%sizepop����ͬʱ������Ǹ��ʴ��ļ����Ŀ
        rand_p = rand(1,1);
        for j=1:archive_num%sizepop����ͬʱ������Ǹ��ʴ��ļ����Ŀ
            if rand_p >= p_fitvalue_ban(j,1) && rand_p < p_fitvalue_ban(j+1,1)
                space_count(j,1) = space_count(j,1)+1;
                break;%ֹͣ��������
            end
        end
    end

%     �����ѡ��
%   ����һ���µ���Ⱥ���Ұ���Ӧ����Ӧֵ�õ�
    leader = zeros(1+6,sizepop);%��pareto������һ����  4����ӵ������
    k = 1;
    for j=1:archive_num
        while space_count(j,1)>0
            leader(:,k) = archive(:,j);%�õ���Ӧ������Ⱥ
            k = k+1;%�±��߶�
            space_count(j,1) = space_count(j,1)-1;%��Ŀ����
        end
    end
    
    disp(leader(1,:));
    %��archive����ĳ�������(��š�������Ӧֵ)֮���ȫ����0 1 5 6 7����Ҫ
    archive((2:1:4),(1:1:2*sizepop)) =  0;%��0
    
    iter = iter+1; 
end

%��������Ӧֵ��ͼ
figure(110);
plot((1:1:ger),aver_fitness_ger(1,:));
hold on;

figure(111);
plot((1:1:ger),aver_fitness_ger(2,:));
hold on;

figure(112);
plot((1:1:ger),aver_fitness_ger(3,:));
hold on;

%��������Ҫ����
save ifa_aver_fitness_ger.mat aver_fitness_ger;%����ƽ��ֵ
save ifa_WSN_fitness.mat new_pop_index_fitness;%���Ľ��
