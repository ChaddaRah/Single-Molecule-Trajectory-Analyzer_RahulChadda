
clear all
close all
%--------------
% 1) Rename the extension of AOI files from .dat to .mat
%
% 2) Add AOIS to empty regions in the Cy3 image and remember how many AOIs were added to estimate bg from emplty feilds
%--------------

% enter the # of AOIs to estimate background
prompt = 'How many AOIs for bg estimation?'
no = input(prompt)


disp ('Select and open Cy3 ROIs')


% Select Imscroll (.mat) file with Cy3 traces
uiopen


% The number of spots 
ntrajectories = aoifits.aoiinfo2(end,6);


% The length of tracjectory
nframes=aoifits.data(end,2);

dataexport=NaN(nframes,(ntrajectories-no)*2);


backgroundCy3= NaN(1,no);
counter=1;

for z=ntrajectories-no+1:ntrajectories
    ts = aoifits.data(z:ntrajectories:end,:);
backgroundCy3(1,counter) = mean(ts(:,8));
counter = counter+1;
end

MeanCy3background = mean(backgroundCy3)

y=[0:ntrajectories-no];

for x=1:ntrajectories-no

ts = aoifits.data(x:ntrajectories:end,:);

dataexport(:,x+y(x))= ts(:,8)-(MeanCy3background); %careful here before adding any number

end




clear aoifits ntrajectories ts x y counter z

disp ('Select and open Cy5 ROIs')

% % Select Imscroll (.mat) file with corresponding Cy5 traces
uiopen


% The number of trajectories in the file
ntrajectories = aoifits.aoiinfo2(end,6);


% The number of frames in the timeseries
nframes=aoifits.data(end,2);

% dataexport=NaN(nframes,ntrajectories*2);
backgroundCy5= NaN(1,no);
counter=1

for z=ntrajectories-no+1:ntrajectories;
    ts = aoifits.data(z:ntrajectories:end,:);
backgroundCy5(1,counter) = mean(ts(:,8));
counter = counter+1;
end

MeanCy5background = mean(backgroundCy5)

y=[1:ntrajectories-no];

for x=1:ntrajectories-no

ts = aoifits.data(x:ntrajectories:end,:);

dataexport(:,x+y(x))= ts(:,8)-MeanCy5background;

end

clear aoifits ts x y

% Write file as tab seperated output 
%writematrix(dataexport,'Cy3-Cy5 traces in vbFRET format (video #1).txt','Delimiter','tab')
%p=dataexport;

%save ('atestMat.mat','p','-ascii')
%save ('atestDat.dat','p','-ascii')

%Write comma limited file for ebFRET calcuation
%dlmwrite('ebFRET INPUT_1to10.dat',dataexport(1:10,1:(ntrajectories-no)*2))

%Preallocate
%Corrected_data = [0 0 0];
% Clean-up Trajectories (remove photobleached regions from traces

counter =1;
for x=1:2:((ntrajectories-no)*2)

datasav = zeros(nframes,8);
    
figure('Renderer', 'painters', 'Position', [50 50 2500 1300]); title(['AOI # ' num2str(x)])
    %Cy3
    subplot(5,1,1); 
    plot(1:nframes,dataexport(:,x),'g','LineWidth',1);
    datacy3=dataexport(:,x); title(['Cy3 AOI # ' num2str(x-1)]) %title ('Cy3')
    %Cy5
    subplot(5,1,2);
    plot(1:nframes,dataexport(:,x+1),'r','LineWidth',1);
    datacy5=dataexport(:,x+1);  title ('Cy5')
    %Cy3/5
    subplot(5,1,3); 
    hold on
    plot(1:nframes,dataexport(:,x),'g','LineWidth',1);
    datacy3=dataexport(:,x);
    plot(1:nframes,dataexport(:,x+1),'r','LineWidth',1);
    datacy5=dataexport(:,x+1); title ('Cy3/5')
    hold off
    %FRET
    subplot(5,1,4); 
    plot(1:nframes ,((datacy5./(datacy5+datacy3))),'m','LineWidth',1);
    ylim([0 1.2]); yticks(0:0.2:1.2); title ('FRET')
    %Total Intensity
    subplot(5,1,5);
    plot(1:nframes,(datacy5+datacy3),'k','LineWidth',1);
    ylim([0 inf]); title ('Total Intensity')
    
    [a,b]= ginput(12);
    a=round(a);
    
    % Retain the portion of trajectory within clicks (rest is set to zero)
%     dataexport((max(a(1),a(2)):end),x)=zeros;
%     dataexport(1:(min(a(1),a(2))),x)=zeros;
%     dataexport((max(a(1),a(2)):end),x+1)=zeros; 
%     dataexport(1:(min(a(1),a(2))),x+1)=zeros;
    
    % IF the selected trace is shorter than 10 frames entire trace is set
    % to zero i.e. not included
%     if length (a(1):a(2)) < 10
%         dataexport(:,x)=zeros;
%         dataexport(:,x+1)=zeros;
%     end 

datasav (1: length(a(1):a(2))',1) = dataexport((a(1):a(2))',x);
datasav (1: length(a(1):a(2))',2) = dataexport((a(1):a(2))',x+1);

datasav (1: length(a(3):a(4))',3) = dataexport((a(3):a(4))',x);
datasav (1: length(a(3):a(4))',4) = dataexport((a(3):a(4))',x+1);

datasav (1: length(a(5):a(6))',5) = dataexport((a(5):a(6))',x);
datasav (1: length(a(5):a(6))',6) = dataexport((a(5):a(6))',x+1);

% datasav (1: length(a(7):a(8))',7) = dataexport((a(7):a(8))',x);
% datasav (1: length(a(7):a(8))',8) = dataexport((a(7):a(8))',x+1);
% 
% datasav (1: length(a(9):a(10))',9) = dataexport((a(9):a(10))',x);
% datasav (1: length(a(9):a(10))',10) = dataexport((a(9):a(10))',x+1);
% 
% datasav (1: length(a(11):a(12))',11) = dataexport((a(11):a(12))',x);
% datasav (1: length(a(11):a(12))',12) = dataexport((a(11):a(12))',x+1);

filename = sprintf('DNA Six (fret roi)_vbFRET INPUT _%d.txt', x);

writematrix(datasav,filename,'Delimiter','tab')
        
    counter = counter+1;
    %w = waitforbuttonpress;
    close
    
    clear a datasav
end

% Remove columns with zeros
%dataexport( :, all(~dataexport,1) ) = [];


% Write file as tab seperated output 
writematrix(datasav,'DNA 6 (fret roi)_vbFRET format_separate traces TAAEST.txt','Delimiter','tab')

% Combine all traces in a single linear array for vbFRET output
% evencolumns = dataexport(:,2:2:end);
% cy5=(evencolumns(:).')';
% cy5 (cy5==0) = [];
% 
% oddcolumns = dataexport(:,1:2:end-1);
% cy3=(oddcolumns(:).')';
% cy3 (cy3==0) = [];
% 
% DonorAcceptor = [cy3 cy5];

%sn = importdata ("F:\RAD52\Rad52 sm data\08-24-2022\Analysis\6. dT62 + ScRad52-delC_1.7 nM (532nm_30 pwr_150 ms + 1 min)\OUTPUT\08-24-2022_1. dT62 + 1.7 nM ScRAd52-DelC_ vbFRET format_separate traces.txt")
%cy3 = sn (:,1);
%cy5 = sn (:,2);

% FRET =(cy5./(cy5 + cy3));
% FRET (FRET<0) = [];

%writematrix(FRET,'08-03-2022_dT80+Rad52_FRET.txt','Delimiter','tab'); % Export FRET traces
% writematrix(DonorAcceptor,'08-03-2022_dT80+Rad52_ vbFRET format_lower_zero_FRET_single concatt.txt','Delimiter','tab'); % Export D/A traces


%% Do you wish to remove noisy FRET points with -ve FRET values??

% A= importdata ('Cy3-Cy5 traces in vbFRET format (composite video)_single concatt.txt');
%   
%     Collect data from different videos
%    %a = importdata('RPA-FAB_A-Cy5-3primeOH_1_ vbFRET format_single concatt.txt');
%    %b = importdata('RPA-FAB_A-Cy5-3primeOH_2_ vbFRET format_single concatt.txt');
%    %c = importdata('RPA-FAB_A-Cy5-3primeOH_3_ vbFRET format_single concatt.txt');
%    %d = importdata('RPA-FAB_A-Cy5-3primeOH_4_ vbFRET format_single concatt.txt');
%    %A = (a,b,c,d);
%
% B =(A(:,2)./(A(:,1) + A(:,2)));
% 
% FRET = horzcat(A,B);
% 
% C = FRET;
% 
% C(C(:, 3)< -0.05, :)= [];
% C(C(:, 3)> 1, :)= [];
% figure
% plot(C(:,3))

% C(:,3)=[];
% 
% writematrix(C,'07-29-2021_SynScRPA-A-Cy5_Cy3-3primeOH_1  (0.35 laser power_100ms)_single concatt (fret outliers removed).txt','Delimiter','tab')

%%


% MeanFRET=NaN(1,ntrajectories-no);
% counter = 1;
% for x=1:2:((ntrajectories-no)*2)
% % Plot intensity vs time
%     figure('Renderer', 'painters', 'Position', [50 50 2500 1300]); title(['AOI # ' num2str(counter)])
%     %Cy3
%     subplot(5,1,1); 
%     plot(1:nframes,dataexport(:,x),'g','LineWidth',1);
%     datacy3=dataexport(:,x); title ('Cy3')
%     %Cy5
%     subplot(5,1,2);
%     plot(1:nframes,dataexport(:,x+1),'r','LineWidth',1);
%     datacy5=dataexport(:,x+1);  title ('Cy5')
%     %Cy3/5
%     subplot(5,1,3); 
%     hold on
%     plot(1:nframes,dataexport(:,x),'g','LineWidth',1);
%     datacy3=dataexport(:,x);
%     plot(1:nframes,dataexport(:,x+1),'r','LineWidth',1);
%     datacy5=dataexport(:,x+1); title ('Cy3/5')
%     hold off
%     %FRET
%     subplot(5,1,4); 
%     plot(1:nframes + length(zeros(50,1)),([(datacy5./(datacy5+datacy3)); zeros(50,1)]),'m','LineWidth',1);
%     ylim([0 1]); yticks(0:0.2:1); title ('FRET')
%     %Total Intensity
%     subplot(5,1,5);
%     plot(1:nframes,(datacy5+datacy3),'k','LineWidth',1);
%     ylim([0 inf]); title ('Total Intensity')
%     
%     FRET = [(datacy5./(datacy5+datacy3)); zeros(50,1)];
%     [a,b]= ginput(2);
%     a=round(a);
%     MeanFRET(counter)=mean(FRET(a(1):a(2)))
%     
%     counter = counter+1;
%     %w = waitforbuttonpress;
%    close
%    
% end
% MeanFRET (MeanFRET<0.1)= []
% MeanFRET (isnan(MeanFRET))= []
% histogram(MeanFRET,20)
% writematrix(MeanFRET','Mean FRET from Traces(video #1).txt','Delimiter','tab')