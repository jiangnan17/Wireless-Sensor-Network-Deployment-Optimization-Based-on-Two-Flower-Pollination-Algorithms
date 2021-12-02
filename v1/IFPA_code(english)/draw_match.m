%%After obtaining the mobile plan through the node planning path algorithm LAPJV algorithm on the fifth page of the paper, redraw the graph.
function draw_match(first_init_sensors,best_indivi,sensor_s)
    global N;
    [match,~] = get_match_value(first_init_sensors,best_indivi,sensor_s);
    %Load the data
    x_pos_init = first_init_sensors(1,:,1);  
    y_pos_init = first_init_sensors(2,:,1);
    x_pos_best = best_indivi(1,:);
    y_pos_best = best_indivi(2,:);
    
    %Used to draw a line that matches the initial position and the final position
    for k=1:N
        plot(x_pos_init(match(k,1)),y_pos_init(match(k,1)),'.','MarkerSize',10);
        plot(x_pos_best(match(k,2)),y_pos_best(match(k,2)),'.','MarkerSize',10);
        
        line([x_pos_init(match(k,1)),x_pos_best(match(k,2))],[y_pos_init(match(k,1)),y_pos_best(match(k,2))]...
               ,'LineWidth',2);
       
        
        %Write the serial number and name of the node 
        %
        if rem(k,4)==1
            %Initial drawing position
            text(x_pos_init(match(k,1))+0.5,y_pos_init(match(k,1)),'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+2,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
            %Optimized drawing position
            text(x_pos_best(match(k,2))+0.5,y_pos_best(match(k,2)),'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+2,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
            
            
        elseif rem(k,4)==2% 
            text(x_pos_init(match(k,1))-2.5,y_pos_init(match(k,1)),'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))-1,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
            
            
            text(x_pos_best(match(k,2))-2.5,y_pos_best(match(k,2)),'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))-1,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        elseif rem(k,4)==3
            %
            text(x_pos_init(match(k,1))- 0.5,y_pos_init(match(k,1))+ 1.5,'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))+0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
        
        
            text(x_pos_best(match(k,2))- 0.5,y_pos_best(match(k,2))+ 1.5,'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))+0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        else% 
            text(x_pos_init(match(k,1))-0.5,y_pos_init(match(k,1))-1.5,'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))-2,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
        
            
            text(x_pos_best(match(k,2))-0.5,y_pos_best(match(k,2))-1.5,'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))-2,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        
        end
        
        hold on;%Add this sentence, otherwise there is one less center
        axis([0,50,0,50]);% 
        set(gca,'xtick',(0:2:50));% 
        set(gca,'ytick',(0:2:50));% 
        axis square;% 
        set(gca,'yminorgrid','on');% 
        set(gca,'xminorgrid','on');
        grid on;% 
        hold on;% 
    end
    

%     Draw a picture of the center point of the grid
%     for i=1:L/1
%         for j=1:W/1
%             plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%             hold on;%Keep the picture
%         end
%     end

        % Draw rhombus
        x1 = [25,35,25,15,25];
        y1 = [35,25,15,25,35];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        hold on;

       % Draw four triangles
        x1 = [0,0,15,0];
        y1 = [50,35,50,50];
        x2 = [0,0,15,0];
        y2 = [15,0,0,15];
        x3 = [35,50,50,35];
        y3 = [50,50,35,50];
        x4 = [35,50,50,35];
        y4 = [0,0,15,0];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        line(x2,y2,'LineWidth',1,'LineStyle','-','color','k');
        line(x3,y3,'LineWidth',1,'LineStyle','-','color','k');
        line(x4,y4,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        fill(x2,y2,[0.5 0.5 0.5]);
        fill(x3,y3,[0.5 0.5 0.5]);
        fill(x4,y4,[0.5 0.5 0.5]);
        hold on;

    

    
end