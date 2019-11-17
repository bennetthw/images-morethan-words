clc, clear all, close all, scrsz = get(groot,'ScreenSize');
% This version is the eleventh version of a script for analyzing cilia or other small objects. 
% From a defined folder of images this script lists samples, creates an
% empty data structure, storing information about the experiment and
% images, with space for ROI coordinates to be stored (in
% "CIP11_Analyze_4manual.m"
% This script can be Run in one step (rather than section by section) and uses dialog boxes to prompt for user input. 
% REQUIREMENTS: Have all the images for
% the chosen experiment in one folder, make a data folder to store output.
% this script analyzes tif images of
% up to 4 channels numbered C0, C1, C2, C3, No stacks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Define input.
% set paths for the location of files to be analyzed.
%https://www.mathworks.com/help/matlab/ref/inputdlg.html
uiwait(msgbox('In the next window, select folder with images you want to analyze.','modal'));
[image_path] = uigetdir('Select folder with images');
experimentpath = ([image_path,'/']); % image path directory

% Generate a temporary sample list for the specified folder
%sample_list_temp=list_samples_everest_temp_projections(experimentpath);
sample_list_temp=list_samples_everest_temp(experimentpath);
for sample_num=1:size(sample_list_temp,1)
    disp(sprintf("\nSample Number: %d, %s", sample_num, char(sample_list_temp(sample_num))));
end
% STEP 2: Define output, and image properties/channel information
% define the folder where data will be stored, and define names.
for sample_num=1:size(sample_list_temp,1)
    indexed_sample_list{sample_num}=['Sample number ',num2str(sample_num),' ',sample_list_temp{sample_num}];
end

uiwait(msgbox('In the next window, select (or create) the folder where data will be stored.','modal'));
[data_path] = uigetdir(experimentpath,'Select folder where data will be stored');
datadir=[data_path,'/'];

Qsample_list=listdlg('PromptString','Select samples to analyze:','ListString',indexed_sample_list,'ListSize',[500 300]);

prompt = {'Enter experimentname:','Enter date:','Channel 1 name', 'Channel 2 name','Channel 3 name','Channel 4 name','Exposure times (List in channel order: DAPI(100ms) FITC(800ms) Cy5(400ms) TexasRed(1000ms) would be 100 800 400 1000','Channel to mask (list mulitple channels to mask a summed image):' 'Channel number for Red in RGB merge, Channel number for Green in RGB merge, Channel number for Blue in RGB merge' 'Left-to-right channel order for mosaics (list 1-4 channel number, separated by spaces):' 'RGB channel assignment for cilia detection image (Put 0 for empty channel):'};
title = 'Input';
dims = [1 100];
definput = {'HWB000_GFP_CEP170_AcTub' datestr(now,'mmddyyyy') 'DAPI', 'FITC' 'Far Red' 'Texas Red' '0 0 0 0' '3' '4 2 3' '4 3 2' '4 0 3'};
answer_step2 = inputdlg(prompt,title,dims,definput);

experimentname=answer_step2{1}; % define a name for this experiment/analysis session
data_initialized_date=char(answer_step2(2)); % also assign a date to all files generated as part of this analysis
%store the list of samples to analyze quantitatively in this session
diary([datadir,experimentname,'_log.txt']); diary on, %Create log file
disp(sprintf("\nRunning CIP12_Initialize_Data_4manual."));

disp(sprintf("\nImages stored in %s", data_path));
disp(sprintf("Data files are written here %s", datadir));
disp(sprintf("File names will contain the experiment name: %s, and date: %s", experimentname,data_initialized_date));
diary off,
disp(sprintf("These samples will be analyzed (Qsample_list): %s", num2str(Qsample_list)));


prompt = {'Channel 1 name','Channel 2 name','Channel 3 name','Channel 4 name','Exposure times (List in channel order: DAPI(100ms) FITC(800ms) Cy5(400ms) TexasRed(1000ms) would be 100 800 400 1000',};
% title = 'Input 3, Image and channel properties';
% dims = [1 100];
% definput = {'DAPI','FITC','Far Red','Texas Red','0 0 0 0'}; % List channel names in the same order as the channels were acquired, usually DAPI (C0), 488(FITC) (C1), 647(cy5)(C2), 568(texas red)(C3)

% answer_step3 = inputdlg(prompt,title,dims,definput);
answer_step3 = {answer_step2{3},answer_step2{4},answer_step2{5},answer_step2{6},answer_step2{7}};
channel_names{1}=answer_step3{1};
channel_names{2}=answer_step3{2};
channel_names{3}=answer_step3{3};
channel_names{4}=answer_step3{4};
exposure_times=str2num(answer_step3{5});

disp(sprintf("\n%s %s", prompt{1},answer_step3{1}));
disp(sprintf("\n%s %s", prompt{2},answer_step3{2}));
disp(sprintf("\n%s %s", prompt{3},answer_step3{3}));
disp(sprintf("\n%s %s", prompt{4},answer_step3{4}));
disp(sprintf("\n%s %s", prompt{5},answer_step3{5}));
% channel_names={'DAPI' 'LAP-hNPHP2 GFP' 'ptg ARL13B' 'CEP170'};
% exposure_times=[160,1200,1000,2000]; % These values have no bearing on the analysis, but it is convenient to have a record of this information.

% indexed_channel_list={['1,  ', channel_names{1}] ['2,  ', channel_names{2}] ['3,  ', channel_names{3}] ['4,  ', channel_names{4}]};
% f=msgbox(indexed_channel_list);
%%
 prompt = {'Channel to mask (list mulitple channels to mask a summed image):' 'Channel number for Red in RGB merge:' 'Channel number for Green in RGB merge:' 'Channel number for Blue in RGB merge:' 'Left-to-right channel order for mosaics (list 1-4 channel number, separated by spaces):' 'RGB channel assignment for cilia detection image (Put 0 for empty channel):'};
answer_step4={answer_step2{8},answer_step2{9},answer_step2{10},answer_step2{11}};

%channel_to_mask=[3]; % choose which channel(s) to use to make a mask
channel_to_mask=str2num(answer_step4{1});
%merge_channel_order=[4 2 3]; % choose R G B channels for merges
merge_channel_order=[str2num(answer_step4{2})]; % choose R G B channels for merges
%channel_order=[4 3 2]; % choose channel order for plots
channel_order=[str2num(answer_step4{3})]; % choose channel order for plots
%ROI_merge_channel_order=[4 0 3]; % choose R G B channels for merges where ROI will be picked, 0 means leave that color blank
ROI_merge_channel_order=[str2num(answer_step4{4})]; % choose R G B channels for merges where ROI will be picked, 0 means leave that color blank

disp(sprintf("\n%s %s", prompt{1},answer_step4{1}));
disp(sprintf("%s %s", prompt{2},num2str(merge_channel_order(1))));
disp(sprintf("%s %s", prompt{3},num2str(merge_channel_order(2))));
disp(sprintf("%s %s", prompt{4},num2str(merge_channel_order(3))));
disp(sprintf("%s %s", prompt{5},answer_step4{3}));
disp(sprintf("Channels will be merged in RGB images for ROI picking in FIJI: %s", num2str(ROI_merge_channel_order)));


% prompt = {'Microns per pixel:' 'Max intenisty:' 'Conversion factor, calculate this as [width of RGB_merge tif in imageJ]/[width of raw tif in pixels] The denominator should be in whatever units are used to report ROI coordinates, the default is inches.'};
% title = 'Input 5, Image properties';
% dims = [1 100];
% definput = {'0.102381' '4095' '1'};
% answer_step5 = inputdlg(prompt,title,dims,definput);
answer_step5 = {'0.102381' '4095' '1'}; % CIP11 exports image in pixels.
% These will most likely be the same for all your experiments.
% microns_per_pixel=0.102381; % Everest 63X oil --> TIFF, can determine with ImageJ image properties, or the log file exported from slidebook
% max_intensity=4095; %4095 for 12 bit images.
% conversion_factor=0.0093; % 0.0093 in/px for Henny's Laptop, imageJ coordinates are in inches, 0.0067 for Apollo Creed, 0.00833 for Apollo Creed considering resolution change
% CIP11 exports image in pixels.conversion factor =1

microns_per_pixel=str2num(answer_step5{1}); % Everest 63X oil --> TIFF, can determine with ImageJ image properties, or the log file exported from slidebook
max_intensity=str2num(answer_step5{2}); %4095 for 12 bit images.
conversion_factor=str2num(answer_step5{3}); % 0.0093 in/px for Henny's Laptop, imageJ coordinates are in inches, 0.0067 for Apollo Creed, 0.00833 for Apollo Creed considering resolution change
disp(sprintf("\n%s %s", prompt{1},answer_step5{1}));
disp(sprintf("%s %s", prompt{2},answer_step5{2}));
disp(sprintf("%s %s", prompt{3},answer_step5{3}));
% Store the permanent version of the  sample list
% the following lines print a report of the information entered above.
% Generate a sample list for the specified folder
if ~exist([datadir,experimentname,'_',data_initialized_date,'_sample_list.mat']);
    list_samples_everest(experimentpath,datadir,experimentname,data_initialized_date);disp(sprintf("The sample list and image reference for this folder have been generated and saved (list_samples_everest.m)"));
else
    disp(sprintf("A pre-existing sample list and image reference for this folder has been found and loaded."));
end
diary on,
load([datadir,experimentname,'_',data_initialized_date,'_sample_list.mat'],'sample_list','image_reference_02'); disp(sprintf("\n%s, Data in date: %s,\n\nSample List:",experimentname,data_initialized_date));
for sample_num=1:size(sample_list,1);
    disp(sprintf("\nSample Number: %d, %s", sample_num, char(sample_list(sample_num))));
end
diary off,
for sample_choice=1:size(sample_list,1);
    % Find the image filenames for the chosen sample
    sample_subset_filenames_everest_nonzero(datadir,experimentname,data_initialized_date,image_reference_02,sample_choice,sample_list);
end
disp(sprintf("\nSample subset filenames listed, sample_subset_filenames_everest_nonzero.m"));

% Step 3. Specify asample from the list, enter channel settings, store image information.
%Prints a report of the channel settings
RGB_name={'Red' 'Green' 'Blue'};descriptive_channel_to_mask=channel_names(channel_to_mask);num_channels_to_plot=size(channel_order,2);
for ch=1:size(channel_order,2);
    descriptive_channel_order(ch)=channel_names(channel_order(ch));
end
num_merge_channel=size(merge_channel_order,2);
for ch=1:size(merge_channel_order,2);
    if merge_channel_order(ch)>0
        descriptive_merge_channel_order(ch)=channel_names(merge_channel_order(ch));
    else
        descriptive_merge_channel_order(ch)={'No channel'};
    end
end
num_ROI_merge_channel=size(ROI_merge_channel_order,2);
for ch=1:size(ROI_merge_channel_order,2);
    if ROI_merge_channel_order(ch)>0
        descriptive_ROI_merge_channel_order(ch)=channel_names(ROI_merge_channel_order(ch));
    else
        descriptive_ROI_merge_channel_order(ch)={'No channel'};
    end
end

uiwait(msgbox('In the next window, adjust the intensity ranges to best display the structure that will be masked. Increase lower bound to reduce background, and decrease upper bound to enhance signal. Hit "OK" to generate a figure.'));
sample_choice=Qsample_list(1);
fieldnum=1;
% DAPI
channel1_min=0.1;
channel1_max=0.5;
% green
channel2_min=0.05;
channel2_max=0.3;
% far red
channel3_min=0.05;
channel3_max=0.3;
% red
channel4_min=0.05;
channel4_max=0.3;

f=msgbox(indexed_sample_list(Qsample_list));
%%
intensities_ok=0;
while ~intensities_ok;
    prompt = {'Sample number to test intensity' 'Field number' 'Channel1 range' 'Channel2 range' 'Channel3 range' 'Channel4 range'};
    title = 'Input 6, intensity determination';
    dims = [1 100];
    definput = {num2str(sample_choice) num2str(fieldnum) [num2str(channel1_min),' ',num2str(channel1_max)] [num2str(channel2_min),' ',num2str(channel2_max)]  [num2str(channel3_min),' ',num2str(channel3_max)]  [num2str(channel4_min),' ',num2str(channel4_max)]};
    answer_step6 = inputdlg(prompt,title,dims,definput); 
    %CHOOSE A SAMPLE (the one you want to determine thresholds with)
    sample_choice=str2num(answer_step6{1}); % Ex. sample_choice=1 chooses the first sample from sample_list
    fieldnum=str2num(answer_step6{2});
    
    % The lines below load the image filenames for the chosen sample
    clear('sample_subset', 'sample_subset_filenames');load([datadir,experimentname,'_',data_initialized_date,'_',num2str(sample_choice),'_',char(sample_list_temp(sample_choice)),'_sample_subset.mat']);diary off;channelnumber=size(fieldnames(sample_subset_filenames),1); % number of channels
    % STEP 3B: Run to test and set Mask Ranges (as many times as needed)
    % Set min, max for each channel to generate mask: Usually trims low (background) signal
    % and saturates high signal.
    % DAPI
    [channel1_ranges]=str2num(answer_step6{3});
    channel1_min=channel1_ranges(1);
    channel1_max=channel1_ranges(2);
    % green
    [channel2_ranges]=str2num(answer_step6{4});
    channel2_min=channel2_ranges(1);
    channel2_max=channel2_ranges(2);
    % far red
    [channel3_ranges]=str2num(answer_step6{5});
    channel3_min=channel3_ranges(1);
    channel3_max=channel3_ranges(2);
    % red
    [channel4_ranges]=str2num(answer_step6{6});
    channel4_min=channel4_ranges(1);
    channel4_max=channel4_ranges(2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mask_ranges=([channel1_min channel1_max; channel2_min channel2_max; channel3_min channel3_max; channel4_min, channel4_max]);
    % Open images for each channel, standardize the intensity range displayed, and make a merge
    Whole_field_subplot_four_channel_CIP08(channelnumber,channel_order,merge_channel_order,experimentpath,sample_subset_filenames,fieldnum,max_intensity,mask_ranges,channel_names,microns_per_pixel,0);
    
    intensity_check=questdlg('Are these intensities OK?','','Yes','No','Yes');
    if strcmp(intensity_check,'No')
        intensities_ok=0;
    elseif strcmp(intensity_check,'Yes')
        intensities_ok=1;
    end
end
diary on,disp(sprintf("\nSamples to analyze: %s\n\nCurrent sample: %d, %s", num2str(Qsample_list),sample_choice, char(sample_list(sample_choice))));disp(sprintf("\n%d Channels", channelnumber));disp(sprintf("%s \n", channel_names{:})); disp(sprintf("Exposure times (ms): "));disp(sprintf("%d\n",exposure_times(:)));disp(sprintf("Max intensity %d", max_intensity));disp(sprintf("microns/pixel %d", microns_per_pixel));disp(sprintf("FIJI to Mat conversion factor: %d", conversion_factor));disp(sprintf("\nMask will be generated from:"));disp(sprintf("%s\n",(descriptive_channel_to_mask{:})));disp(sprintf("Up to 4 channels can be plotted, %d have been specified. \nImages will be plotted in this order:",num_channels_to_plot));disp(sprintf("%s\n",(descriptive_channel_order{:})));disp(sprintf("3 channels can be merged, %d have been specified. \nChannels will be merged in this order when plotting figures:",num_merge_channel));disp(sprintf("%s: %s",(RGB_name{1}),(descriptive_merge_channel_order{1})));disp(sprintf("%s: %s",(RGB_name{2}),(descriptive_merge_channel_order{2})));disp(sprintf("%s: %s",(RGB_name{3}),(descriptive_merge_channel_order{3})));disp(sprintf("\n3 channels can be merged, %d have been specified. \nChannels will be merged in this order when saving images for ROI selection:",num_ROI_merge_channel));disp(sprintf("%s: %s",(RGB_name{1}),(descriptive_ROI_merge_channel_order{1})));disp(sprintf("%s: %s",(RGB_name{2}),(descriptive_ROI_merge_channel_order{2})));disp(sprintf("%s: %s",(RGB_name{3}),(descriptive_ROI_merge_channel_order{3})));
disp(sprintf("\nChannel 1, %s range for object detection and masking: %s %s",channel_names{1},num2str(channel1_min), num2str(channel1_max)));
disp(sprintf("\nChannel 2, %s range for object detection and masking: %s %s",channel_names{2},num2str(channel2_min), num2str(channel2_max)));
disp(sprintf("\nChannel 3, %s range for object detection and masking: %s %s",channel_names{3},num2str(channel3_min), num2str(channel3_max)));
disp(sprintf("\nChannel 4, %s range for object detection and masking: %s %s",channel_names{4},num2str(channel4_min), num2str(channel4_max)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 5. Save a images and Initialize data
    export_check=questdlg('Ready to export images?','','Yes','No','Yes');
    if strcmp(export_check,'No')
        export_ok=0;
    elseif strcmp(export_check,'Yes')
        export_ok=1;
    end
    
if export_ok
    %%
    close all
    f=msgbox('Exporting images...');
    stackimagedir=[datadir,'Stacks/']; mkdir(stackimagedir);
    for sample_choice=Qsample_list
        % Load the image filenames for the chosen sample
        clear('sample_subset', 'sample_subset_filenames');
        load([datadir,experimentname,'_',data_initialized_date,'_',num2str(sample_choice),'_',char(sample_list(sample_choice)),'_sample_subset.mat']);
        diary off
        
        imdir=[datadir,'ManAdjMerges/']; 
        %For separate folders for each sample uncomment the line below
        %imdir=[datadir,'ManAdjMerges_Sample_',num2str(sample_choice),'/'];
        if ~isdir(imdir)
            mkdir(imdir);
        end
        stackimagedir=[datadir,'Stacks/','Sample_',num2str(sample_choice),'/']; mkdir(stackimagedir);
        for fieldnum=1:size(sample_subset_filenames,2) %Do it for all fields
            [imname]=adaptive_cilia_merge_tif_CIP12(ROI_merge_channel_order,experimentpath,sample_subset_filenames,fieldnum,max_intensity,mask_ranges,imdir,stackimagedir,sample_subset);
            data(sample_choice).FIJI_tif(fieldnum).name=[imdir,imname,'.tif'];
            data(sample_choice).FIJI_csv(fieldnum).name=[imname,'.csv'];
            data(sample_choice).path=experimentpath;
            data(sample_choice).date=data_initialized_date;
            data(sample_choice).Qsample_list=Qsample_list;
            data(sample_choice).sample(fieldnum).image=sample_subset_filenames(fieldnum);
            data(sample_choice).sample(fieldnum).image.cilia_coordinates=zeros(1,6);
            data(sample_choice).shifts(fieldnum).X_shift_channel2=0;
            data(sample_choice).shifts(fieldnum).Y_shift_channel2=0;
            data(sample_choice).shifts(fieldnum).X_shift_channel4=0;
            data(sample_choice).shifts(fieldnum).Y_shift_channel4=0;
            data(sample_choice).mask_ranges=mask_ranges;
            data(sample_choice).channel_names=channel_names;
            data(sample_choice).exposure_times=exposure_times;
            data(sample_choice).conversion_factor=conversion_factor;
            data(sample_choice).microns_per_pixel=microns_per_pixel;
            data(sample_choice).max_intensity=max_intensity;
            data(sample_choice).channel_to_mask.num=channel_to_mask;
            data(sample_choice).channel_to_mask.name=descriptive_channel_to_mask;
            data(sample_choice).channel_order.num=channel_order;
            data(sample_choice).channel_order.name=descriptive_channel_order;
            data(sample_choice).merge_channel_order.num=merge_channel_order;
            data(sample_choice).merge_channel_order.name=descriptive_merge_channel_order;
            data(sample_choice).ROI_merge_channel_order.num=ROI_merge_channel_order;
            data(sample_choice).ROI_merge_channel_order.name=descriptive_ROI_merge_channel_order;
        end
    end
    save([datadir,experimentname,'_',data_initialized_date,'_data.mat'],'data');
    close(f);
    diary on, disp(sprintf("Images were exported using adaptive_cilia_merge_tif_CIP12, and data structure was initialized. Open images in FIJI and save coordinates in csv files")); diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%