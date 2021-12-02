%%Draw the minimum spanning tree connected to the network
%The smallest spanning tree algorithm is the Kruskal algorithm
function draw_MST(x_pos,y_pos,sersors_r,adjacencyMatrix, adjacencyMatrix_dis)
    global N;
    angle=0:pi/100:2*pi;%
    for k=1:N
          %The nodes of different radii are processed with different colors
        if sersors_r(1,k) == 7
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','r');% 
            hold on;
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','r');% 
            hold on;
        elseif sersors_r(1,k) == 6
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','g');% 
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','g');% 
            hold on;
        else
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','b');% 
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','b');% 
            hold on;
        end
        
        if rem(k,4)==1
            text(x_pos(1,k)+0.5,y_pos(1,k),'N','FontSize',15);% 
            text(x_pos(1,k)+2,y_pos(1,k)-0.5,num2str(k),'FontSize',10);% 
        elseif rem(k,4)==2% 
            text(x_pos(1,k)-2.5,y_pos(1,k),'N','FontSize',15);%
            text(x_pos(1,k)-1,y_pos(1,k)-0.5,num2str(k),'FontSize',10);% 
        elseif rem(k,4)==3
            % 
            text(x_pos(1,k)- 0.5,y_pos(1,k)+ 1.5,'N','FontSize',15);% 
            text(x_pos(1,k)+1,y_pos(1,k)+0.5,num2str(k),'FontSize',10);% 
        else% 
            text(x_pos(1,k)-0.5,y_pos(1,k)-1.5,'N','FontSize',15);% 
            text(x_pos(1,k)+1,y_pos(1,k)-2,num2str(k),'FontSize',10);% 
        end
        hold on;% 
        
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

       %Draw four triangles
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

    
    %Make the final drawing minimum spanning tree
%     load adjacencyMatrix_dis.mat;
%     load adjacencyMatrix.mat;
%     disp('Weights');
%     disp(adjacencyMatrix_dis);
%     disp('Whether there is an edge');
%     disp(adjacencyMatrix);
%     pause(10);
    [weight_sum, span_tree] = kruskal(adjacencyMatrix, adjacencyMatrix_dis);
    [row,~] = size(span_tree);%Get the horizontal and vertical coordinates of span_tree
%     disp('row');
%     disp(row);
%     disp('colu');
%     disp(colu);
%     disp('span_tree');
%     disp(span_tree);

    %Pay attention to the writing of line
    for i=1:row
        line([x_pos(1,span_tree(i,1)),x_pos(1,span_tree(i,2))],[y_pos(1,span_tree(i,1)),y_pos(1,span_tree(i,2))]...
           ,'LineWidth',2);
    end
    %disp(['The total weight of the minimum spanning tree£º',num2str(weight_sum)]);
    hold on;
end