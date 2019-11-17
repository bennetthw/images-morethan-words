% This script loads the data file. Imports ROI coordinates. Stores mosaic images, makes masks
% and measurements.

clc, clear all, close all
uiwait(msgbox('In the next window, select the data file for the experiment you want to analyze.','modal'));
[data_filename, data_pathname] = uigetfile;
load([data_pathname,data_filename]);

% check paths
[update_path]=check_paths_data(data_pathname,data);
if ~isempty(update_path)
    for i=1:size(data,2)
        if ~isempty(data(i).path)
            data(i).path=update_path;
        end
    end
end

datadir=data_pathname;
for i=1:size(data,2)
    if ~isempty(data(i).path)
        experimentpath=data(i).path;
        experimentname1=strsplit(data_filename,['_',data(i).date,'_data.mat']);
        experimentname=experimentname1{1};
        date=data(i).date;
        Qsample_list=data(i).Qsample_list;
    end
end
%
load([datadir,experimentname,'_',date,'_sample_list.mat'],'sample_list','image_reference_02');
load('CIP12_Analyze_Indexed_Data_Names.mat'); data_number_guide=msgbox(indexed_data_names);
diary([datadir,experimentname,'_log.txt']); diary on,
disp(sprintf("\nRunning CIP12_Analyze."));


%% INPUT
use_classifier=0;
mask_parameters.channel_to_mask=data(Qsample_list(1)).channel_to_mask.num;
mask_parameters.threshMan=0.3;
mask_parameters.blurring=5;
mask_parameters.mask_min=40; %minimum object size in pixels
mask_parameters.dilation=1;
mask_parameters.background_distance=6;
mask_parameters.background_dilation=3;
mask_parameters.mask_ranges=data(Qsample_list(1)).mask_ranges;

channel_to_mask=mask_parameters.channel_to_mask;
threshMan=mask_parameters.threshMan;
blurring=mask_parameters.blurring;
dilation=mask_parameters.dilation;
background_distance=mask_parameters.background_distance;
background_dilation=mask_parameters.background_dilation;
blurring=mask_parameters.blurring;
mask_ranges=mask_parameters.mask_ranges;

replicate_number=1;
tilesize=[13 2.5];  %tile width and tile height in microns
sort_data_choice=15; %this is the data parameter that will be used for sorting, choose 0 to not sort.
channel_toweight=4;

%Prompt to edit parameters
%
prompt = {'Enter channel(s) to mask (Listing multiple channels here creates a mask of the sum of those channels. Additional masks can be added later.)','Enter mask threshold (0-1, usually ~0.3, meaning that positive signal is 30% above background)','Bottom weighted channel (ROIs will be rotated so the end with the highest signal in this channel is the bottom)','Tile size ([width height] of thumbnail in microns)','Replicate number (to save independent files for measurements made with different parameters)'};
title = 'Input for masking';
dims = [1 100];
definput = {num2str(channel_to_mask),num2str(threshMan),num2str(channel_toweight),num2str(tilesize),num2str(replicate_number)};
answer_stepM = inputdlg(prompt,title,dims,definput);
channel_to_mask=str2num(answer_stepM{1});
for i=Qsample_list
    data(i).channel_to_mask.num=channel_to_mask;
    data(i).channel_to_mask.name=data(i).channel_names(channel_to_mask);
end
save([data_pathname,data_filename],'data');
    
threshMan=str2num(answer_stepM{2});
channel_toweight=str2num(answer_stepM{3});
tilesize=round((str2num(answer_stepM{4}))./data(Qsample_list(1)).microns_per_pixel);
replicate_number=str2num(answer_stepM{5});

% %%%Full prompt,
% prompt = {'Use data from classifier? (0 for no, 1 for yes)', 'Enter channels to mask','Enter mask theshold (0-1, usually ~0.3)','Blurring factor','Mask dilation (pixels)','Distance between object mask and background mask (pixels)','Thickness of background mask (pixels)','Replicate number (to save independent files for measurements made with different parameters','Tile size ([width height] of thumbnail in pixels)','Measurement to sort by: (15 is cilia length)','Bottom weighted channel'};
% title = 'Input for masking';
% dims = [1 100];
% definput = {num2str(use_classifier),num2str(channel_to_mask),num2str(threshMan),num2str(blurring),num2str(dilation),num2str(background_distance),num2str(background_dilation),num2str(replicate_number),num2str(tilesize),num2str(sort_data_choice),num2str(channel_toweight)};
% answer_stepM = inputdlg(prompt,title,dims,definput);
% use_classifier=str2num(answer_stepM{1});
% channel_to_mask=str2num(answer_stepM{2});
% threshMan=str2num(answer_stepM{3});
% blurring=str2num(answer_stepM{4});
% dilation=str2num(answer_stepM{5});
% background_distance=str2num(answer_stepM{6});
% background_dilation=str2num(answer_stepM{7});
% replicate_number=str2num(answer_stepM{8});
% tilesize=str2num(answer_stepM{9});
% sort_data_choice=str2num(answer_stepM{10});
% channel_toweight=str2num(answer_stepM{11});


disp(sprintf("Masking and Measuring parameters \n Analysis replicate number %d",replicate_number));
disp(sprintf("Channel(s) to mask %s",num2str(channel_to_mask)));
disp(sprintf("Masking threshold %d", threshMan)); 
disp(sprintf("Blurring %d", blurring)); 
disp(sprintf("Mask dilation %d", dilation)); 
disp(sprintf("Background mask distance %d", background_distance)); 
disp(sprintf("Background mask thickness %d", background_dilation));
disp(sprintf("Data sorted by measurement No. %d", sort_data_choice));
disp(sprintf("Tile width: %d", tilesize(2))); 
disp(sprintf("Tile height: %d", tilesize(1))); 
disp(sprintf("Channel enriched at ciliary base (channel to weight) %d",channel_toweight));

analyzedatetime=datestr(now,'mmddyyyyTHHMM');
% % Dialog box: Make a folder where training images and xml files will be
% saved, if exporting directly.
% uiwait(msgbox('In the next window, select (or create) the folder where images for the classifier will be stored. THIS MUST BE IN PUBLIC FOLDER.','modal'));
% [trainclassifierimg_path] = uigetdir('Volumes/public/','Select folder where data will be stored');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if use_classifier;
    % Import cilia coordinates from Nandita's csv
    disp(sprintf("Analyzing ROIs from classifier."));
    uiwait(msgbox('In the next window, select the folder containing the csvs.','modal'));
    [csv_pathname] = uigetdir;
    csvdir=[csv_pathname,'/'];
    Import_Data_from_csv_CIP10(datadir,csvdir,experimentname,date,Qsample_list);
    f=msgbox('Object ROI coordinates imported from csv file.');
    % Stores mosaics using Nandita's coordinates.
    f=msgbox('Storing ROIs...');
    Mosaic_Data_Saves_CIP10_4classifier(datadir,experimentname,date,date,Qsample_list,tilesize);
    close(f);
    f=msgbox('Validating ROIs...');
    Mosaic_Data_Saves_CIP10_classifierROIvalidator(datadir,experimentname,date,date,Qsample_list);
    f=msgbox('ROI images saved.');
    close(f);
    f=msgbox('Masking_and_measuring_01_CIP10.m Making masks and measurements...');
    Masking_and_measuring_01_CIP11_imported_data(datadir,experimentname,date,Qsample_list,mask_parameters,replicate_number);
    close(f);
    f=msgbox('Measurements made.');
    close(f);
else %prompts to import coordinates form FIJI csvs, in which case a folder containing these files must be specified. Otherwise expects the coordinates to have already been stored in the data variable
    disp(sprintf("Analyzing ROIs from manual determination."));
    uiwait(msgbox('In the next window, select the folder containing the csvs exported from FIJI.','modal'));
    [csv_pathname] = uigetdir;
    if csv_pathname ~=0 %Import the coordinates
        csvdir=[csv_pathname,'/'];
        Import_Data_from_csv_CIP11(datadir,csvdir,experimentname,date,Qsample_list);
        f=msgbox('Object ROI coordinates imported from csv file.');
        close(f);
        f=msgbox('Storing ROIs...');
        % Stores mosaic data structures for all samples
        Mosaic_Data_Saves_CIP10(datadir,experimentname,date,date,Qsample_list,tilesize);
        close(f);
        f=msgbox('Validating ROIs...');
        Mosaic_Data_Saves_CIP10_AnnotateROI(datadir,experimentname,date,date,Qsample_list,tilesize);
        close(f);
        f=msgbox('ROI images saved.');
        close(f);
        f=msgbox('Masking_and_measuring_01_CIP12_imported_data.m Making masks and measurements...');
        Masking_and_measuring_01_CIP12_imported_data(datadir,experimentname,date,Qsample_list,mask_parameters,replicate_number);
        close(f);
        f=msgbox('Measurements made.');
        close(f);
        
    else %coordinates were not imported from csv files, the masking and measuring functions will use the data file open in the workspace, rather than loading one.
        f=msgbox('Storing ROIs...');
        % Stores mosaic data structures for all samples
        Mosaic_Data_Saves_CIP11(data,datadir,experimentname,date,date,Qsample_list,tilesize);
        close(f);
        f=msgbox('Validating ROIs...');
        Mosaic_Data_Saves_CIP11_AnnotateROI(data,sample_list,datadir,experimentname,date,date,Qsample_list,tilesize);
        close(f);
        f=msgbox('ROI images saved.');
        close(f);
        f=msgbox('Masking_and_measuring_01_CIP12.m Making masks and measurements...');
        Masking_and_measuring_01_CIP12(data,datadir,experimentname,date,Qsample_list,mask_parameters,replicate_number);
        close(f);
        f=msgbox('Measurements made.');
        close(f);
    end
end
%
% Sort by measurment
Sort_only_CIP10(datadir,experimentname, date, Qsample_list,replicate_number, sort_data_choice);

% Combine all measurements into one file.
subset=['Rep',num2str(replicate_number)];
Combine_only_CIP10(datadir,experimentname,Qsample_list, subset, date, date);

% Load data
subset=['Rep',num2str(replicate_number)];
load([datadir,experimentname,'_',date,'_sample_list.mat'],'sample_list','image_reference_02');
load([datadir,experimentname,'_',subset,'_',date,'_combined_data.mat'],'combined_sorted_data');

exportCSDcsv(combined_sorted_data,Qsample_list)

% Rotate the cilia
mask_ranges=combined_sorted_data(Qsample_list(1)).mosaic_data_sorted.mask_ranges;         
min_int1=mask_ranges(1,1); max_int1=mask_ranges(1,2); %DAPI, channel1
min_int2=mask_ranges(2,1); max_int2=mask_ranges(2,2); %Green channel, channel2
min_int3=mask_ranges(3,1); max_int3=mask_ranges(3,2); %Blue channel, channel3
min_int4=mask_ranges(4,1); max_int4=mask_ranges(4,2); %Red channel, channel4
samples=Qsample_list; midwayfactor=0.5; determine_bottom_automatically=1; %choose 1, to use the channel to weight
for sample_choice=Qsample_list
    [RotatedCSD1,rotatedate]=CIP11_Rotate_Cilia_multi_mask(combined_sorted_data, data_pathname, channel_toweight, midwayfactor, determine_bottom_automatically, sample_choice, tilesize, min_int2, max_int2, min_int3, max_int3, min_int4, max_int4);
    RotatedCSD1(sample_choice).mosaic_data_sorted.Qsample_list=Qsample_list;
    RotatedCSD(sample_choice).mosaic_data_sorted=RotatedCSD1(sample_choice).mosaic_data_sorted;
    save([datadir,experimentname,'_',rotatedate,'_rotated_data_SampleNumber_',num2str(sample_choice),'.mat'],'-v7.3','RotatedCSD1');
    disp(sprintf("Saved Rotated Data, %s", [datadir,experimentname,'_',rotatedate,'_rotated_data_SampleNumber_',num2str(sample_choice),'.mat']));
end
disp(sprintf("COMPLETED, CIP12 ANALYZE."));
%%
% % Open this script to correct the channel alignment.
% CIP10_correct_translation.m
diary off