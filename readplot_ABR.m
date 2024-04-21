%% Program to read all exported ABRs, create an age-specific grand averaged trace with STD  
% Reads data files from given folders, stores header info, plots time domain waveform.
% Made to work with CFTS format from new electrophys software/system (TDT)
% written by AP and LSA on 4/15/24 based on previous code 

%Tweaked for ABR's by Luz Andrino -- April 2024 
%trying to loop through all files in a folder to average across then plot
clear;
clc;


%% Readme: this section is important and controls analysis/directory. 

in_data_dir= 'C:\Users\ASP97\University of Pittsburgh\SHRS Partha Lab - Documents\Data\Animal_Data\Analysis\AgingSeries_Gerbils\dAM Aging Exports\75wk\ABR\3kHz\'; % path to binary files (input)
data_out_dir= 'C:\Users\ASP97\University of Pittsburgh\SHRS Partha Lab - Documents\Data\Animal_Data\Analysis\AgingSeries_Gerbils\dAM Aging Data\75wk\ABR\'; % where output ABR values will be saved 
% two files will be saved: a mat file, and a figure (png)
fig_out_dir= data_out_dir; % saving the output file in the same folder 

do_save_data= 1; % whether to save data 
do_save_fig= 1; % whether to save figure 

allfiles= dir([in_data_dir 'GER*_E1']); % assuming that all 

fs = 24414.0625; %TDT sampling frequency 
filttype = 1; % 1 for filtering
lp = 3000; hp = 80; order=4;  %Check if this is OK

%create a matrix for grand average data (time domain) 
allfiles_tdatsum = [];
allfiles_tdatdiff = [];

%%
car_freq_Hz= 3e3;
num_of_files= length(allfiles);
all_names= {allfiles.name}';

for fileVar=1:length(allfiles)
    fileVar
        cur_fStruct= allfiles(fileVar);
        fname= [cur_fStruct.folder filesep cur_fStruct.name];
        [level, tdattemp, WF] = read_cfts_abr_data_AP(fname);
        arr_size = size(tdattemp);
        tdatsumtemp = tdattemp(1:(arr_size(1)/2),:);     %First half of the response exported from CFTS is the sum of two polarities
        tdatdifftemp = tdattemp((arr_size(1)/2)+1:arr_size(1),:); %Second half of the response is the difference of two polarities
        
        %filter the data 
        if filttype == 1
          
        [bh,ah]=butter(order,hp/(fs*0.5),'high');
        [bl,al] = butter(order,lp/(fs*0.5),'low');
        i = 0;
    for i = 1:arr_size(2)
        thfiltsum(:,i) = filtfilt(bh,ah,tdatsumtemp(:,i));
        tlfiltsum(:,i) = filtfilt(bl, al, thfiltsum(:,i));
     end
    tdatsumtemp = tlfiltsum;
        end 
        
        
        
        %Get data ready for grand average
                   
        allfiles_tdatsum = [allfiles_tdatsum, tdatsumtemp(:,1:9)];
         allfiles_tdatdiff = [allfiles_tdatdiff, tdatdifftemp(:,1:9)];% all 9 levels are appended next to each other for all animals here 
%         %call every 9th column when constructing GA later. This is a shitty way to do it. change later. 
%             allfiles_tdat = [allfiles_tdat, tdattemp(:,1:9)];
        tdatsumtemp = [];
        tdatdifftemp = [];

end 

%% create GA for each sound level
for j = 1:9
ABR_GA_sum(:,j) = mean(allfiles_tdatsum(:,j:9:size(allfiles_tdatsum,2)),2);
ABR_GA_diff(:,j) = mean(allfiles_tdatdiff(:,j:9:size(allfiles_tdatdiff,2)),2);
end 

%% split file for first 20ms 



%% apply filter 



%% plots

figure 
for i = 1:9 
    subplot(9,1,i)
    plot(ABR_GA_sum(:,i))
    xlim ([0 250])
end 


%%

% for k = 1:9
% for j = k:9:size(allfiles_tdat,2)
%     
   

%Apply grand average 
% grand_avg = mean(allfiles_tdat,2);


%%

% 
% %Loop
% for i = 1:length(file_list)
%     fn = fullfile("D:\Exports_G1921\ABR\500Hz\10step", file_list(i).name);
% 
% %fn = 'C:\Codigo\GER1322 ABR Exports Pre\GER13221-1_E1';
% fs = 24414.0625; %TDT sampling frequency 
% 
% filttype = 1; % 1 for filtering
% 
% zp = 1; % multiple for zero padding
% 
% freq = 500; % enter frequency of AM 
% CarrierFreq = 1000; % enter carrier frequency
% 
% [level, tdattemp, WF] = read_cfts_abr_data_AP(fn);
% 
% %Get data ready for grand average
% allfiles_tdatsum = [allfiles_tdatsum, tdattemp(:,1)];
% 
% % plot(tvec, tdattemp(:, 1), 'b');
% % hold on;
% 
% %Apply grand average 
% grand_avg = mean(allfiles_tdatsum,2);
% 
% %Plot please??
% figure; 
% tvec = (1: size(grand_avg,1))/fs;
% 
% 
% 
% end
% 
% 


%%
% 
% tdat(:,1) = tdattemp(:,1);   % forD:\ABR Exports\G1222\ABR\Pre\GER12221_1_E2 future batch processing
% arr_size = size(tdat);
% 
% % Plot Averaged Time traces 
% 
% tvec=1/fs:1/fs:arr_size(1)/fs;%Time domain vector at sampling frequency
% tdat =tdat';
% figure
% subplot(3,1,1)
%  plot(tvec,tdat,'b')
%  title('Timetrace from ABR CFTS file', 'Fontsize',11)
%  
%  tdatsum = tdat(1, 1:(arr_size(1)/2));     %First half of the response exported from CFTS is the sum of two polarities
%  tdatdiff = tdat(1, (arr_size(1)/2)+1:arr_size(1)); %Second half of the response is the difference of two polarities, not needed in this case
%  
%  tvecsum = tvec(1, 1:(arr_size(1)/2));
%  subplot(3,1,2)
%  plot(tvecsum,tdatsum,'b')
%  subplot(3,1,3)
%  plot(tvecsum,tdatdiff,'b')
 

%%

 