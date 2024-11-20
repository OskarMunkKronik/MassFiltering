
%% Run Massfiltering pipeline
%% Inputs: 
% Data: MS1 and MS2 data of an LCxLC-HRMS measurement as .mat files
% mzroi_aug: m/z-features extracted with ROI 
% Rt: retention times
% Modulationtime
%
%% Parameters:
% 1.Peak Picking:
% threshold_peak = Minimum Peak Height --> defined in options_MF.txt
% threhsold_mass = Allowed mass deviation --> defined in options_MF.txt
% 
% 2.Mass Filtering:
% threshold_sim_1 = Similarity threshold when MF is used alone --> defined in options_MF.txt
% threshold_sim_2 = Similarity threhsold when MF as pre-processing --> defined in options_MF.txt
% xwindow = retention time window in D2 for peak extraction --> defined in options_MF.txt
% ywindow = retention time window in D1 for peak extraction --> defined in options_MF.txt
%
% 3. Deconvolution (MCR / SIT)
% MaxIter = maximum number of iterations --> defined in options_SIT.txt
% ConvCrit = Convergence criterion --> defined in options_SIT.txt
% Constr_1 = non-negativity in elution profiles (1 = on / 0 = off) --> defined in options_SIT.txt
% Constr_2 = non-negativity in mass spectra (1 = on / 0 = off) --> defined in options_SIT.txt
% Init = Initialization (0 = random) --> defined in options_SIT.txt
%
% 4. Parameters not yet included in init file
% cutoff = mz-cutoff for discarding mz-values > M+H + cutoff

%% Requirements: 
% Matlab Tensor toolbox (https://www.tensortoolbox.org/)
% Data base look up is currently only possible with local .MSP files

display('select project folder')
path_to_project_folder = uigetdir();
addpath(path_to_project_folder)
cd(path_to_project_folder)
addpath(genpath('Code'))
load('Data\Data_all.mat');
%% Creating the reshaped three-way data set from the unfolded data set
MFC = MassFiltering_class(); % initializing the class
MFC = MFC.prepare(Data_new,mzroi_aug,Rt_new,Modulationtime); % Reshaping the data
%% Initializing the DataAnalysis class
DAC = DataAnalysis_class(); % initializing the class
DAC = DAC.init(MFC); % defining project repository, loading option files and suspect list

%% Running the peak picking based on the exact masses in the suspect list
screening = 1; % 1 = automatically select the highest peak as M+1, 0 = manually select the M+1 peak
DAC = DAC.check_peaks(screening); % execute peak picking

%% Run Mass Filtering
filtered = 0; % 0 = MF without additional MCR
cutoff = 50; % masses that are larger than M+1+cutoff Da are removed
DAC = DAC.RunMassFiltering(filtered,cutoff); % Run Mass Filtering
filtered = 1; % 1 = MF with additional MCR
cutoff = 50; % masses that are larger than M+1+cutoff Da are removed
DAC = DAC.RunMassFiltering(filtered,cutoff); % Run Mass Filtering

%% Run MCR
filtered = 0; % 0 = Run MCR on original data, no pre-filtering
DAC = DAC.get_fac(filtered); % get the rank for determining the number of components in MCR
DAC = DAC.run_mcr(filtered); % run MCR

%%
filtered = 1; % 1 = Run MCR on pre-filtered data
model_type = "MCR"; % Defines whether "MCR" or "SIT" should be used as model
DAC = DAC.get_fac(filtered); % get the rank for determining the number of components in MCR
DAC = DAC.run_mcr(filtered); % run MCR

%% Extract the resulting mass spectra and save them to the repository

DAC = DAC.get_results();
%% Extract data base spectra from Mass Bank (https://github.com/MassBank/MassBank-data/releases)

% !!! this function is not optimized and currently uses the suspect name to
% search in Mass Bank for entries with the same name !!! --> slow and prone
% to fail
databaselookup2()

%% Creating three files summarizing the similarity of the reference spectrum with the extracted spectra and the number of shared fragments
DAC = DAC.statistical_analysis();
