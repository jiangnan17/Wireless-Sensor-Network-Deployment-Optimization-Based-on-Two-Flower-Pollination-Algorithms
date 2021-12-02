%%��Բ
function draw_MST(x_pos,y_pos,sersors_r,adjacencyMatrix, adjacencyMatrix_dis)
    global N;
    angle=0:pi/100:2*pi;%�Ƕ�
    for k=1:N
        %�����²�ͬ�뾶��ɫ��ͬ
        if sersors_r(1,k) == 7
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','r');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            hold on;
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','r');%��Բ
            hold on;
        elseif sersors_r(1,k) == 6
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','g');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','g');%��Բ
            hold on;
        else
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','b');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','b');%��Բ
            hold on;
        end
        
        if rem(k,4)==1
            text(x_pos(1,k)+0.5,y_pos(1,k),'N','FontSize',15);%���б��
            text(x_pos(1,k)+2,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%���б��
        elseif rem(k,4)==2%�ұ�
            text(x_pos(1,k)-2.5,y_pos(1,k),'N','FontSize',15);%���б��
            text(x_pos(1,k)-1,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%���б��
        elseif rem(k,4)==3
            %�ϱ�
            text(x_pos(1,k)- 0.5,y_pos(1,k)+ 1.5,'N','FontSize',15);%���б��
            text(x_pos(1,k)+1,y_pos(1,k)+0.5,num2str(k),'FontSize',10);%���б��
        else%�±�
            text(x_pos(1,k)-0.5,y_pos(1,k)-1.5,'N','FontSize',15);%���б��
            text(x_pos(1,k)+1,y_pos(1,k)-2,num2str(k),'FontSize',10);%���б��
        end
        hold on;%�������  ��Ȼ��һ��Բ��
        
        axis([0,50,0,50]);%����������С�����ֵ  
        set(gca,'xtick',(0:2:50));%����x���경��Ϊ1
        set(gca,'ytick',(0:2:50));%����y���경��Ϊ1   ���и����� ����������̫������
        axis square;%ʹ���ݱ���Ϊ1  ����Բ����Բ  ��Ȼ����Բ
        set(gca,'yminorgrid','on');%��С����
        set(gca,'xminorgrid','on');
        grid on;%����
        hold on;%��ͼһֱ��������
    end
    
    % ���������ĵ��ͼ������
%     for i=1:L/1
%         for j=1:W/1
%             plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%             hold on;%��ͼһֱ��������
%         end
%     end

        %������1
        x1 = [15,35,35,15,15];
        y1 = [40,40,35,35,40];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        hold on;

        %������2
        x2 = [22,28,28,22,22];
        y2 = [35,35,15,15,35];
        line(x2,y2,'LineWidth',1,'LineStyle','-','color','k');
        fill(x2,y2,[0.5 0.5 0.5]);
        hold on;

        %������3
        x3 = [15,35,35,15,15];
        y3 = [15,15,10,10,15];
        line(x3,y3,'LineWidth',1,'LineStyle','-','color','k');
        fill(x3,y3,[0.5 0.5 0.5]);
        hold on;

        %���ĸ�������
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

    
    %%�������Ļ���С������
%     load adjacencyMatrix_dis.mat;
%     load adjacencyMatrix.mat;
%     disp('Ȩ��');
%     disp(adjacencyMatrix_dis);
%     disp('�Ƿ���ڱ�');
%     disp(adjacencyMatrix);
%     pause(10);
    [weight_sum, span_tree] = kruskal(adjacencyMatrix, adjacencyMatrix_dis);
    [row,~] = size(span_tree);%�õ�span_tree�ĺ�������
%     disp('row');
%     disp(row);
%     disp('colu');
%     disp(colu);
%     disp('span_tree');
%     disp(span_tree);
    %ע��line��д��
    for i=1:row
        line([x_pos(1,span_tree(i,1)),x_pos(1,span_tree(i,2))],[y_pos(1,span_tree(i,1)),y_pos(1,span_tree(i,2))]...
           ,'LineWidth',2);
    end
    %disp(['��С����������Ȩֵ��',num2str(weight_sum)]);
    hold on;
end