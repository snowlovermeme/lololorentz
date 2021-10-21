for U=5500:500:8000
N=5;
interval="5000";  % 放电间隔，单位毫秒

Tca=[];Tcm=[];
Tia=[];Tim=[];
h = waitbar(0,'仿真进程');
model.param.set('U0', num2str(U));
for i=1:N
    s=['第' num2str(i) '次放电:初始温升阶段'];
    waitbar(i/N,h,s)
    
    
    % 求解器设置
if i==1
    model.study("std3").feature("time").set("useinitsol", false);
else
    model.study("std3").feature("time").set("useinitsol", true);
    model.study("std3").feature("time").set("initmethod", "sol");
    model.study("std3").feature("time").set("initstudy", "std4");
    model.study("std3").feature("time").set("solnum", "last");
    model.study("std3").feature("time").set("solvertype", "solnum");
    model.study("std3").feature("time").set("solnumhide", "off");
    model.study("std3").feature("time").set("initstudyhide", "on");
    model.study("std3").feature("time").set("initsolhide", "on");
end
% 设置放电间隔
model.study("std4").feature("time").set("tlist", "range(0,5,"+interval+")");
%  运行
model.study("std3").run();
    s=['第' num2str(i) '次放电：冷却阶段'];
    waitbar(i/N,h,s)
model.study("std4").run();

    s=['第' num2str(i) '次放电：处理数据阶段'];
    waitbar(i/N,h,s)

% 后处理，导出数据
% 每一匝铜的平均温度
model.result().table("tbl4").clearTableData();
model.result().numerical("av1").setResult();
for j=1:12  
      
    model.result().numerical("av1").selection().set((j-1)*2+6);
    model.result().numerical("av1").set("table", "tbl4");
    model.result().numerical("av1").appendResult();  
end
tab1 = mphtable(model,"tbl4").data(:,2:end);
Tca=[Tca;tab1];

% 计算每一匝绝缘的平均温度
model.result().table("tbl4").clearTableData();
model.result().numerical("av1").setResult();
for j=1:13    
    model.result().numerical("av1").selection().set((j-1)*2+5);
    model.result().numerical("av1").set("table", "tbl4");
    model.result().numerical("av1").appendResult(); 
end
tab2 = mphtable(model,'tbl4').data(:,2:end);
Tia=[Tia;tab2];

% 计算放电结束后每一匝的最大值
model.result().table("tbl5").clearTableData();
model.result().numerical("max1").setResult();

for j=1:12
    model.result().numerical("max1").selection().set((j-1)*2+6);
    model.result().numerical("max1").set("dataseries", "maximum");
    model.result().numerical("max1").set("data", "dset4");
    model.result().numerical("max1").set("table", "tbl5")
    model.result().numerical("max1").appendResult();
end
tab3 = mphtable(model,"tbl5").data;
Tcm=[Tcm;tab3];
%=============================================================%
model.result().table("tbl5").clearTableData();
model.result().numerical("max1").setResult();
h1=waitbar(0,"最大值计算");
for j=1:13
    waitbar(j/13,h1);
model.result().numerical("max1").selection().set((j-1)*2+5);
model.result().numerical("max1").set("dataseries", "maximum");
model.result().numerical("max1").set("data", "dset4");
    model.result().numerical("max1").set("table", "tbl5");
    model.result().numerical("max1").appendResult();
end
tb4=mphtable(model,"tbl5").data;
Tim=[Tim;tab3];

close(h1)
tab4 = mphtable(model,"tbl5").data;

% 导出gif
model.result().export("anim1").set("giffilename", "J:\simulationModel\comsol\model\多次放电的温升\数据存储\gif"+num2str(i)+".gif");

model.result().export("anim1").run();
end
% [m,n]=size(Tcm);
% time=1:1:m;
% turn=1:13;
% surf(time,turn,Tcm')
% shading interp
% figure(1)
% hold on
% for i=1:13
%     plot(time,Tca(:,i))
%     legend(num2str(i))
% end

xlswrite("Tca-"+num2str(U)+".xls", Tca);
xlswrite("Tia-"+num2str(U)+".xls", Tia);
xlswrite("Tcm-"+num2str(U)+".xls", Tcm);
xlswrite("Tim-"+num2str(U)+".xls", Tim);

end
close(h)
