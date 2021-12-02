function createhggroup(Parent1, XData1, YData1, XData2, YData2, XData3, YData3, XData4, YData4, XData5, YData5, YData6, YData7, YData8, YData9, YData10, XData6, YData11, XData7, YData12, XData8, YData13, XData9, YData14, XData10, YData15, YData16, YData17, YData18, YData19, YData20, XData11, YData21, XData12, YData22, XData13, YData23, XData14, YData24, XData15, YData25, XData16, YData26, XData17, YData27, XData18, YData28, XData19, YData29, XData20, YData30, XData21, YData31, XData22, YData32, XData23, YData33, XData24, YData34, XData25, YData35)
%CREATEHGGROUP(PARENT1, XDATA1, YDATA1, XDATA2, YDATA2, XDATA3, YDATA3, XDATA4, YDATA4, XDATA5, YDATA5, YDATA6, YDATA7, YDATA8, YDATA9, YDATA10, XDATA6, YDATA11, XDATA7, YDATA12, XDATA8, YDATA13, XDATA9, YDATA14, XDATA10, YDATA15, YDATA16, YDATA17, YDATA18, YDATA19, YDATA20, XDATA11, YDATA21, XDATA12, YDATA22, XDATA13, YDATA23, XDATA14, YDATA24, XDATA15, YDATA25, XDATA16, YDATA26, XDATA17, YDATA27, XDATA18, YDATA28, XDATA19, YDATA29, XDATA20, YDATA30, XDATA21, YDATA31, XDATA22, YDATA32, XDATA23, YDATA33, XDATA24, YDATA34, XDATA25, YDATA35)
%  PARENT1:  hggroup parent
%  XDATA1:  line xdata
%  YDATA1:  line ydata
%  XDATA2:  line xdata
%  YDATA2:  line ydata
%  XDATA3:  line xdata
%  YDATA3:  line ydata
%  XDATA4:  line xdata
%  YDATA4:  line ydata
%  XDATA5:  line xdata
%  YDATA5:  line ydata
%  YDATA6:  line ydata
%  YDATA7:  line ydata
%  YDATA8:  line ydata
%  YDATA9:  line ydata
%  YDATA10:  line ydata
%  XDATA6:  line xdata
%  YDATA11:  line ydata
%  XDATA7:  line xdata
%  YDATA12:  line ydata
%  XDATA8:  line xdata
%  YDATA13:  line ydata
%  XDATA9:  line xdata
%  YDATA14:  line ydata
%  XDATA10:  line xdata
%  YDATA15:  line ydata
%  YDATA16:  line ydata
%  YDATA17:  line ydata
%  YDATA18:  line ydata
%  YDATA19:  line ydata
%  YDATA20:  line ydata
%  XDATA11:  line xdata
%  YDATA21:  line ydata
%  XDATA12:  line xdata
%  YDATA22:  line ydata
%  XDATA13:  line xdata
%  YDATA23:  line ydata
%  XDATA14:  line xdata
%  YDATA24:  line ydata
%  XDATA15:  line xdata
%  YDATA25:  line ydata
%  XDATA16:  line xdata
%  YDATA26:  line ydata
%  XDATA17:  line xdata
%  YDATA27:  line ydata
%  XDATA18:  line xdata
%  YDATA28:  line ydata
%  XDATA19:  line xdata
%  YDATA29:  line ydata
%  XDATA20:  line xdata
%  YDATA30:  line ydata
%  XDATA21:  line xdata
%  YDATA31:  line ydata
%  XDATA22:  line xdata
%  YDATA32:  line ydata
%  XDATA23:  line xdata
%  YDATA33:  line ydata
%  XDATA24:  line xdata
%  YDATA34:  line ydata
%  XDATA25:  line xdata
%  YDATA35:  line ydata

%  由 MATLAB 于 04-Jul-2019 09:24:32 自动生成

% 创建 hggroup
hggroup1 = hggroup('Parent',Parent1);

% 创建 line
line(XData1,YData1,'Parent',hggroup1,'Tag','Upper Whisker','LineStyle','--');

% 创建 line
line(XData2,YData2,'Parent',hggroup1,'Tag','Upper Whisker','LineStyle','--');

% 创建 line
line(XData3,YData3,'Parent',hggroup1,'Tag','Upper Whisker','LineStyle','--');

% 创建 line
line(XData4,YData4,'Parent',hggroup1,'Tag','Upper Whisker','LineStyle','--');

% 创建 line
line(XData5,YData5,'Parent',hggroup1,'Tag','Upper Whisker','LineStyle','--');

% 创建 line
line(XData1,YData6,'Parent',hggroup1,'Tag','Lower Whisker','LineStyle','--');

% 创建 line
line(XData2,YData7,'Parent',hggroup1,'Tag','Lower Whisker','LineStyle','--');

% 创建 line
line(XData3,YData8,'Parent',hggroup1,'Tag','Lower Whisker','LineStyle','--');

% 创建 line
line(XData4,YData9,'Parent',hggroup1,'Tag','Lower Whisker','LineStyle','--');

% 创建 line
line(XData5,YData10,'Parent',hggroup1,'Tag','Lower Whisker',...
    'LineStyle','--');

% 创建 line
line(XData6,YData11,'Parent',hggroup1,'Tag','Upper Adjacent Value');

% 创建 line
line(XData7,YData12,'Parent',hggroup1,'Tag','Upper Adjacent Value');

% 创建 line
line(XData8,YData13,'Parent',hggroup1,'Tag','Upper Adjacent Value');

% 创建 line
line(XData9,YData14,'Parent',hggroup1,'Tag','Upper Adjacent Value');

% 创建 line
line(XData10,YData15,'Parent',hggroup1,'Tag','Upper Adjacent Value');

% 创建 line
line(XData6,YData16,'Parent',hggroup1,'Tag','Lower Adjacent Value');

% 创建 line
line(XData7,YData17,'Parent',hggroup1,'Tag','Lower Adjacent Value');

% 创建 line
line(XData8,YData18,'Parent',hggroup1,'Tag','Lower Adjacent Value');

% 创建 line
line(XData9,YData19,'Parent',hggroup1,'Tag','Lower Adjacent Value');

% 创建 line
line(XData10,YData20,'Parent',hggroup1,'Tag','Lower Adjacent Value');

% 创建 line
line(XData11,YData21,'Parent',hggroup1,'Tag','Box','Color',[0 0 1]);

% 创建 line
line(XData12,YData22,'Parent',hggroup1,'Tag','Box','Color',[0 0 1]);

% 创建 line
line(XData13,YData23,'Parent',hggroup1,'Tag','Box','Color',[0 0 1]);

% 创建 line
line(XData14,YData24,'Parent',hggroup1,'Tag','Box','Color',[0 0 1]);

% 创建 line
line(XData15,YData25,'Parent',hggroup1,'Tag','Box','Color',[0 0 1]);

% 创建 line
line(XData16,YData26,'Parent',hggroup1,'Tag','Median','Color',[1 0 0]);

% 创建 line
line(XData17,YData27,'Parent',hggroup1,'Tag','Median','Color',[1 0 0]);

% 创建 line
line(XData18,YData28,'Parent',hggroup1,'Tag','Median','Color',[1 0 0]);

% 创建 line
line(XData19,YData29,'Parent',hggroup1,'Tag','Median','Color',[1 0 0]);

% 创建 line
line(XData20,YData30,'Parent',hggroup1,'Tag','Median','Color',[1 0 0]);

% 创建 line
line(XData21,YData31,'Parent',hggroup1,'Tag','Outliers',...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','+',...
    'LineStyle','none');

% 创建 line
line(XData22,YData32,'Parent',hggroup1,'Tag','Outliers',...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','+',...
    'LineStyle','none');

% 创建 line
line(XData23,YData33,'Parent',hggroup1,'Tag','Outliers',...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','+',...
    'LineStyle','none');

% 创建 line
line(XData24,YData34,'Parent',hggroup1,'Tag','Outliers',...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','+',...
    'LineStyle','none');

% 创建 line
line(XData25,YData35,'Parent',hggroup1,'Tag','Outliers',...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','+',...
    'LineStyle','none');

