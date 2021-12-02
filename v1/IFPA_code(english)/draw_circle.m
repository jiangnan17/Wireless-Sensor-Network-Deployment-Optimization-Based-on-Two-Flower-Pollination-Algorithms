%This function is used to draw circles and obstacles, etc., for displaying images
function draw_circle(x_pos,y_pos,sersors_r)
    global N;
    angle=0:pi/100:2*pi;%angle
    for k=1:N
        %The nodes of different radii are processed with different colors
        if sersors_r(1,k) == 7
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','r');% Draw the center of the circle, first draw the center of the circle, then draw the circle, otherwise there will be one less
            hold on;
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','r');%Draw a circle
            hold on;
        elseif sersors_r(1,k) == 6
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','g');
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','g');
            hold on;
        else
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','b');
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','b');
            hold on;
        end
        
        %Used to display the serial number and name of the node, up, down, left and right are divided into four directions
        if rem(k,4)==1
            text(x_pos(1,k)+0.5,y_pos(1,k),'N','FontSize',15);% To mark
            text(x_pos(1,k)+2,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%
        elseif rem(k,4)==2%right
            text(x_pos(1,k)-2.5,y_pos(1,k),'N','FontSize',15);%
            text(x_pos(1,k)-1,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%
        elseif rem(k,4)==3
            %up
            text(x_pos(1,k)- 0.5,y_pos(1,k)+ 1.5,'N','FontSize',15);%进行标记
            text(x_pos(1,k)+1,y_pos(1,k)+0.5,num2str(k),'FontSize',10);%进行标记
        else%down
            text(x_pos(1,k)-0.5,y_pos(1,k)-1.5,'N','FontSize',15);%进行标记
            text(x_pos(1,k)+1,y_pos(1,k)-2,num2str(k),'FontSize',10);%进行标记
        end
        hold on;%Add this sentence, otherwise there is one less center
        
        axis([0,50,0,50]);%Abscissa and max coordinates  
        set(gca,'xtick',(0:2:50));%Set the x coordinate step to 2
        set(gca,'ytick',(0:2:50));%
        axis square;%Make the horizontal to vertical ratio 1 so that the circle looks like a circle, otherwise it looks like an ellipse
        set(gca,'yminorgrid','on');%Minimal grid
        set(gca,'xminorgrid','on');
        grid on;%Gridding
        hold on;%Keep the picture
    end
    
%     Draw a picture of the center point of the grid
%     for i=1:L/1
%         for j=1:W/1
%             plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%             hold on;%Keep the picture
%         end
%     end

    %Draw rhombus
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
end