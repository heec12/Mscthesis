fileID = fopen('nf0m05ppp.txt','r');
formatSpec = '%f %f';
sizeV = [2 Inf];
V = fscanf(fileID,formatSpec,sizeV);
fclose(fileID);
V = V';
l = length(V);

t = [0:l-1];
t = t';
V = [t V];
Vleft = V(:,2);
Vright = V(:,3);

one = ones(1,l);

Fs = 1000;
L = length(Vright);

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',20,'HalfPowerFrequency2',500, ...
               'DesignMethod','butter','SampleRate',Fs);
fvtool(d,'Fs',Fs)
vvright = filtfilt(d,Vright);
vvleft = filtfilt(d,Vleft);
t = [1:L];

plot(t,Vright,t,Vleft'color',[0,0,0]+0.6,'Linewidth',2)
hold on;
plot(t,vvright,t,vvleft,'color',[0.75,0,0],'Linewidth',1)
ylim([-4 4]);
%%xticks([20 40 60 80 100 120 140 160 180 200])
%%xticklabels({'100','200','300','400','500','600','700','800','900','1000'})
ylabel('Mean horizontal velocity')
xlabel('Time (kyr)')
%title('kinematic M=0.8')
grid

%newV = [vvleft vvright]; 
%fileIDD = fopen('./filt_fbm08.txt','w');
%fprintf(fileIDD,'%f %f\n',newV);
%fclose(fileIDD);