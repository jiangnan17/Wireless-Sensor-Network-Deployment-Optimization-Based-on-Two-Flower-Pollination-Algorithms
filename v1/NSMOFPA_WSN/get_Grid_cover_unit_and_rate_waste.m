%%得到覆盖率
function [cover_rate,waste_rate] = get_Grid_cover_unit_and_rate_waste(sensor_mat1,sersors_r,per_sersons_num,per_sersons_radius_type)
    %浪费的比例 = （理想覆盖 - 真实覆盖 - 重叠覆盖）/ 理想覆盖
    global M;
    global N;
    global Grid_cen_x_and_y;
    global L;
    global W;
    
    %最笨的办法 分别写出来
    %进行数据的处理  通过截取
    fiftynum = (1:1:50);
    caseone_row1 = [fiftynum(1:1:15),fiftynum(36:1:50)];
    caseone_row2 = [fiftynum(1:1:14),fiftynum(37:1:50)];
    caseone_row3 = [fiftynum(1:1:13),fiftynum(38:1:50)];
    caseone_row4 = [fiftynum(1:1:12),fiftynum(39:1:50)];
    caseone_row5 = [fiftynum(1:1:11),fiftynum(40:1:50)];
    caseone_row6 = [fiftynum(1:1:10),fiftynum(41:1:50)];
    caseone_row7 = [fiftynum(1:1:9),fiftynum(42:1:50)];
    caseone_row8 = [fiftynum(1:1:8),fiftynum(43:1:50)];
    caseone_row9 = [fiftynum(1:1:7),fiftynum(44:1:50)];
    caseone_row10 = [fiftynum(1:1:6),fiftynum(45:1:50)];
    caseone_row11 = [fiftynum(1:1:5),fiftynum(46:1:50)];
    caseone_row12 = [fiftynum(1:1:4),fiftynum(47:1:50)];
    caseone_row13 = [fiftynum(1:1:3),fiftynum(48:1:50)];
    caseone_row14 = [fiftynum(1:1:2),fiftynum(49:1:50)];
    caseone_row15 = [fiftynum(1:1:1),fiftynum(50:1:50)];



    caseone_row36 = [fiftynum(1:1:1),fiftynum(50:1:50)];
    caseone_row37 = [fiftynum(1:1:2),fiftynum(49:1:50)];
    caseone_row38 = [fiftynum(1:1:3),fiftynum(48:1:50)];
    caseone_row39 = [fiftynum(1:1:4),fiftynum(47:1:50)];
    caseone_row40 = [fiftynum(1:1:5),fiftynum(46:1:50)];
    caseone_row41 = [fiftynum(1:1:6),fiftynum(45:1:50)];
    caseone_row42 = [fiftynum(1:1:7),fiftynum(44:1:50)];
    caseone_row43 = [fiftynum(1:1:8),fiftynum(43:1:50)];
    caseone_row44 = [fiftynum(1:1:9),fiftynum(42:1:50)];
    caseone_row45 = [fiftynum(1:1:10),fiftynum(41:1:50)];
    caseone_row46 = [fiftynum(1:1:11),fiftynum(40:1:50)];
    caseone_row47 = [fiftynum(1:1:12),fiftynum(39:1:50)];
    caseone_row48 = [fiftynum(1:1:13),fiftynum(38:1:50)];
    caseone_row49 = [fiftynum(1:1:14),fiftynum(37:1:50)];
    caseone_row50 = [fiftynum(1:1:15),fiftynum(36:1:50)];
    Grid_cover_unit = zeros(L,W);%1*M个用于保存2500个网格的联合概率
    Grid_cover_bool = zeros(L,W,N);%每个网格中心监测的概率  只能是0或者
    %这里也得考虑  涉及到求重叠区域
    for i=1:L
		for k=1:W
            
            %笨的方法处理四个三角形
            if i == 1 && (any(caseone_row1==k))==1
                continue;
            end
            if i == 2 && (any(caseone_row2==k))==1
                continue;
            end

            if i == 3 && (any(caseone_row3==k))==1
                continue;
            end
            if i == 4 && (any(caseone_row4==k))==1
                continue;
            end
            if i == 5 && (any(caseone_row5==k))==1
                continue;
            end
            if i == 6 && (any(caseone_row6==k))==1
                continue;
            end
            if i == 7 && (any(caseone_row7==k))==1
                continue;
            end
            if i == 8 && (any(caseone_row8==k))==1
                continue;
            end
            if i == 9 && (any(caseone_row9==k))==1
                continue;
            end
            if i == 10 && (any(caseone_row10==k))==1
                continue;
            end
            if i == 11 && (any(caseone_row11==k))==1
                continue;
            end
            if i == 12 && (any(caseone_row12==k))==1
                continue;
            end
            if i == 13 && (any(caseone_row13==k))==1
                continue;
            end
            if i == 14 && (any(caseone_row14==k))==1
                continue;
            end
            if i == 15 && (any(caseone_row15==k))==1
                continue;
            end


            if i == 36 && (any(caseone_row36==k))==1
                continue;
            end
            if i == 37 && (any(caseone_row37==k))==1
                continue;
            end
            if i == 38 && (any(caseone_row38==k))==1
                continue;
            end
            if i == 39 && (any(caseone_row39==k))==1
                continue;
            end
            if i == 40 && (any(caseone_row40==k))==1
                continue;
            end
            if i == 41 && (any(caseone_row41==k))==1
                continue;
            end
            if i == 42 && (any(caseone_row42==k))==1
                continue;
            end
            if i == 43 && (any(caseone_row43==k))==1
                continue;
            end
            if i == 44 && (any(caseone_row44==k))==1
                continue;
            end
            if i == 45 && (any(caseone_row45==k))==1
                continue;
            end
            if i == 46 && (any(caseone_row46==k))==1
                continue;
            end

            if i == 47 && (any(caseone_row47==k))==1
                continue;
            end
            if i == 48 && (any(caseone_row48==k))==1
                continue;
            end
            if i == 49 && (any(caseone_row49==k))==1
                continue;
            end

            if i == 50 && (any(caseone_row50==k))==1
                continue;
            end

            %矩形障碍物上的判断
            if (i>35&&i<=40)&&(k>15&&k<=35)%处理掉矩形的
                continue;
            end
            %矩形障碍物中的判断
            if (i>15&&i<=35)&&(k>22&&k<=28)%处理掉矩形的
                continue;
            end
            %矩形障碍物下的判断
            if (i>10&&i<=15)&&(k>15&&k<=35)%处理掉矩形的
                continue;
            end
            
			for j=1:N
				if ((Grid_cen_x_and_y(i,k,1)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(i,k,2)-sensor_mat1(2,j))^2)...
						<=sersors_r(1,j)^2
					Grid_cover_bool(i,k,j) = 1;%代表第i网格中心点可以被第j个传感器节点覆盖
				end
			end
		end
    end
    
    
    %%计算联合分布概率
    for i=1:L
		for k=1:W
            
            %笨的方法处理四个三角形
            if i == 1 && (any(caseone_row1==k))==1
                continue;
            end
            if i == 2 && (any(caseone_row2==k))==1
                continue;
            end

            if i == 3 && (any(caseone_row3==k))==1
                continue;
            end
            if i == 4 && (any(caseone_row4==k))==1
                continue;
            end
            if i == 5 && (any(caseone_row5==k))==1
                continue;
            end
            if i == 6 && (any(caseone_row6==k))==1
                continue;
            end
            if i == 7 && (any(caseone_row7==k))==1
                continue;
            end
            if i == 8 && (any(caseone_row8==k))==1
                continue;
            end
            if i == 9 && (any(caseone_row9==k))==1
                continue;
            end
            if i == 10 && (any(caseone_row10==k))==1
                continue;
            end
            if i == 11 && (any(caseone_row11==k))==1
                continue;
            end
            if i == 12 && (any(caseone_row12==k))==1
                continue;
            end
            if i == 13 && (any(caseone_row13==k))==1
                continue;
            end
            if i == 14 && (any(caseone_row14==k))==1
                continue;
            end
            if i == 15 && (any(caseone_row15==k))==1
                continue;
            end


            if i == 36 && (any(caseone_row36==k))==1
                continue;
            end
            if i == 37 && (any(caseone_row37==k))==1
                continue;
            end
            if i == 38 && (any(caseone_row38==k))==1
                continue;
            end
            if i == 39 && (any(caseone_row39==k))==1
                continue;
            end
            if i == 40 && (any(caseone_row40==k))==1
                continue;
            end
            if i == 41 && (any(caseone_row41==k))==1
                continue;
            end
            if i == 42 && (any(caseone_row42==k))==1
                continue;
            end
            if i == 43 && (any(caseone_row43==k))==1
                continue;
            end
            if i == 44 && (any(caseone_row44==k))==1
                continue;
            end
            if i == 45 && (any(caseone_row45==k))==1
                continue;
            end
            if i == 46 && (any(caseone_row46==k))==1
                continue;
            end

            if i == 47 && (any(caseone_row47==k))==1
                continue;
            end
            if i == 48 && (any(caseone_row48==k))==1
                continue;
            end
            if i == 49 && (any(caseone_row49==k))==1
                continue;
            end

            if i == 50 && (any(caseone_row50==k))==1
                continue;
            end

            %矩形障碍物上的判断
            if (i>35&&i<=40)&&(k>15&&k<=35)%处理掉矩形的
                continue;
            end
            %矩形障碍物中的判断
            if (i>15&&i<=35)&&(k>22&&k<=28)%处理掉矩形的
                continue;
            end
            %矩形障碍物下的判断
            if (i>10&&i<=15)&&(k>15&&k<=35)%处理掉矩形的
                continue;
            end
            
			p = 1;%设置为1
			for j=1:N
				p = p*(1-Grid_cover_bool(i,k,j));%根据公式计算
			end
			Grid_cover_unit(i,k) = 1-p;%联合覆盖概率
		end
    end

   %计算覆盖率
    p_sum = 0;%累加所有网格的覆盖概率
    for i=1:L
        for j=1:W
            p_sum = p_sum + Grid_cover_unit(i,j);
        end
    end
    
    
    cover_rate = p_sum*(L*W/M)/((L*W)- 480 - 320);%覆盖率%留意这里，其实这里折中是480    并不是450 根据点来计算为512
    
    
    %真实的覆盖面积
    area_real = p_sum;
    %disp(['真实面积：',num2str(area_real)]);
    
    %统计重叠区域
    area_over = 0;%重叠面积
    Grid_cover_over_matri = zeros(L,W);
    for i=1:L
        for j=1:W
            Grid_cover_over_matri(i,j) = sum(Grid_cover_bool(i,j,:));%统计个监测点被覆盖多少次
        end
    end
    %计算每个监测点被重复覆盖的次数
    for i=1:L
        for j=1:W
            if Grid_cover_over_matri(i,j) >= 2
                area_over  = area_over + Grid_cover_over_matri(i,j) - 1;%重复 n-1次
            end
        end
    end
    
    %disp(['重叠面积：',num2str(area_over)]);
    
    
    %利用网格点来计算，不然有可能出现负数 因为网格点有些误差
    r_max_area = 158;
    r_mid_area = 117;
    r_min_area = 81;
    %存放理想的覆盖面积
    area_ideal = r_max_area*per_sersons_num(1,1) + r_mid_area*per_sersons_num(1,2) + r_min_area*per_sersons_num(1,3);
    
    

    %disp(['理想面积：',num2str(area_ideal)]);
    waste_rate = (area_ideal - area_real - area_over)/area_ideal;
end