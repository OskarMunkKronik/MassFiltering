
classdef MassFiltering_class
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Reshaped_Data
        rt_axis
        mz_temp
        modulationlength
        modulationtime
        TIC_RT2_Data
        all_intervals
        ind_slices_to_cut
    end

    methods

        function obj = prepare(obj,Data,mzroi,Rt,modulationtime,ind)
            %I
            %Data,mzroi,Rt,modulationtime,rt2_window
%             obj = varargin{1};
%             Data = varargin{2};
%             mzroi = varargin{4};
%             Rt = varargin{5};
%             modulationtime = varargin{6};
            if exist('ind','var')
                if length(ind) == 2
                    rt2_window = ind;
                    clear ind;
                end
            end

            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            data_temp{1} = Data{1,1};
            data_temp{2} = Data{1,2};

            drop_masses{1} = sum(data_temp{1},2)==0;
            data_temp{1}   = data_temp{1}(~drop_masses{1},:);
            data_temp_tic_ms{1} = sum(data_temp{1},1);

            drop_masses{2} = sum(data_temp{2},2)==0;
            data_temp{2}   = data_temp{2}(~drop_masses{2},:);
            data_temp_tic_ms{2} = sum(data_temp{2},1);

            % define all properties of object
            obj.mz_temp{1} = mzroi(~drop_masses{1});
            obj.TIC_RT2_Data{1} = sum(data_temp{1},1);

            obj.mz_temp{2} = mzroi(~drop_masses{2});
            obj.TIC_RT2_Data{2} = sum(data_temp{2},1);

            obj.Reshaped_Data = [];
            obj.rt_axis = Rt;
            obj.modulationtime = modulationtime;
            obj.modulationlength = [];
            obj.all_intervals.rt1 = [];
            obj.all_intervals.rt2 = [];
            obj.ind_slices_to_cut = [];
            % interactively select a search window if rt2_window has not
            % been defined
            if ~exist('rt2_window','var')
                close all
                f = figure;
                plot(obj.TIC_RT2_Data{1})
                while ishandle(f)
                    %prompt for x
                    xlimInput = input('Enter new X limits [xmin xmax]: ');
                    %Check if the input is valid
                    if numel(xlimInput) == 2
                        % set the new limits
                        xlim(xlimInput)
                    else
                        display('Invalid input.')
                    end
                end
                rt2_window = xlimInput;
            end
            % creating the reshaped data set
            for kk = 1:2
                data_new = {};
                rt_new = {};
                masses = {};
                i = 1;
                ii = 1;

                if kk == 1

                    while length(data_temp_tic_ms{1}) > i+rt2_window(2)
                        ind_slices_to_cut(ii) = i;
                        ind_min = find(data_temp_tic_ms{1}(i+rt2_window(1):i+rt2_window(2)) == min(data_temp_tic_ms{1}(i+rt2_window(1):i+rt2_window(2))));
                        data_new{ii} = data_temp_tic_ms{1}(i:i+rt2_window(1)+ind_min);
                        rt_new{ii}   = obj.rt_axis(i:i+rt2_window(1)+ind_min);
                        masses{ii}   = data_temp{1}(:,i:i+rt2_window(1)+ind_min);
                        ii = ii+1;
                        i = i+rt2_window(1)+ind_min+1;
                    end
                    ind_slices_to_cut(ii) = i;
                else
                    ind = ind_slices_to_cut;
                    for k = 1:length(ind)-1
                        data_new{ii} = data_temp_tic_ms{2}(ind(ii):ind(ii+1)-1);
                        rt_new{ii}   = obj.rt_axis(ind(ii):ind(ii+1)-1);
                        masses{ii}   = data_temp{2}(:,ind(ii):ind(ii+1)-1);
                        ii = ii+1;
                        i = ind(end);
                    end

                end

                data_new{ii} = data_temp_tic_ms{kk}(i:end);
                rt_new{ii} =    obj.rt_axis(i:end);
                masses{ii}   = data_temp{kk}(:,i:end);
                for iii = 1:numel(data_new)
                    data_sizes(iii) = numel(data_new{iii});
                end
                ind_max = max(data_sizes);
                data_temp{kk} = [];
                for iii = 1:numel(data_new)
                    if ind_max > numel(data_new{iii})
                        dif = ind_max-numel(data_new{iii});
                        obj.modulationlength(iii) = numel(data_new{iii});
                        fix_value = data_new{iii}(end);
                        fix_mass  = masses{iii}(:,end);
                        data_new{iii} = [data_new{iii} repmat(fix_value,1,dif)];
                        masses{iii} = [masses{iii} reshape(repmat(fix_mass,1,dif),length(fix_mass),dif)];
                    end
                    if iii == 1
                        ind1 = [];
                        ind2 = [];
                        ind3 = [];
                        vals = [];
                    end
                    [ind] = find(masses{iii} ~= 0);
                    [indtemp1, indtemp2] = ind2sub(size(masses{iii}),ind);
                    ind1(length(ind1)+1:length(ind1)+length(indtemp1)) = indtemp1;
                    vals(length(ind2)+1:length(ind2)+length(indtemp2)) = masses{iii}(masses{iii} ~= 0);
                    ind2(length(ind2)+1:length(ind2)+length(indtemp2)) = indtemp2;

                    ind3(length(ind3)+1:length(ind3)+length(indtemp1)) = repmat(iii,1,length(indtemp1));


                end

                tensorsize = [ind_max,numel(data_new),length(obj.mz_temp{kk})];
                subscripts = [ind2; ind3; ind1]';

                obj.Reshaped_Data{kk} = sptensor(subscripts,vals',tensorsize);
                try
                    obj.ind_slices_to_cut = ind_slices_to_cut;
                catch
                    obj.ind_slices_to_cut = 0;
                end
            end
        end


        function obj = run_massfiltering(obj,mz1,axis_mz1,mz2,axis_mz2,options)

            %% Inputs

            % suspect list containing suspect names and M+H masses in separate columns
            % mz1                     = sparse_tensor containing the reshaped MS-1 data.
            % axis_mz1                = feature list of mz-fragments for MS-1
            % mz2                     = sparse_tensor containing the reshaped MS-2 data.
            % axis_mz2                = feature list of mz-fragments for MS-2

            % options.threshold_peak   = threshold for picking peaks in LC1-LC2 plane for the target M+H+. Every peak with intensity below that value will not be selected as possible precursor Ion.
            % options.threshold_masses = value between 0 and 1 that defines the selectivity for selecting masses for the reconstruction of the compound spectrum
            % options.directory        = directory for saving the results
            % options.ID               = either "sample" or "standard" or "blank"
       

            %set default if not defined
            if isempty(options.threshold_peak)
                threshold_peak = 1e04;
            else
                threshold_peak = options.threshold_peak;
            end

            if isempty(options.threshold_masses)
                threshold_masses = 0.9;
            else
                threshold_masses = options.threshold_masses;
            end

            if isempty(options.directory)
                path = pwd;
            else
                path = options.directory;
                cd(path)
            end


            % Load suspect list
            [filename,filepath] = uigetfile('*.*','select a spreadsheet file.');
            disp(['Selected file:' fullfile(filepath,filename)]);

            %read selected file
            suspect_list = readtable(fullfile(filepath,filename))

            %define column names for M/Z-Information
            column_name_mz = input('define column name for M/Z-Information: '); % uncomment this
            eval(strcat('suspect_info_masses = suspect_list.',column_name_mz))

            %define column names for analyte names
            column_name_id = input('define column name for analyte names: '); % uncomment this
            eval(strcat('suspect_info_id = suspect_list.',column_name_id))
            sanitizeString = @(str) regexprep(str,{'-',',',';','[',']','(',')',' ','+'},'');
            suspect_info_id = cellfun(sanitizeString,suspect_info_id,'UniformOutput',false);
            for i = 2:max(size(suspect_info_id))

                mz_target = suspect_info_masses(i);

                [~,ind_mz1] = min(abs(axis_mz1-mz_target));
                [~,ind_mz2] = min(abs(axis_mz2-mz_target));

                if length(suspect_info_id{i}) > 20
                    suspectname = suspect_info_id{i}(1:20);
                else
                    suspectname = suspect_info_id{i};
                end

                foldername = convertStringsToChars(strcat(suspectname,'_',options.ID));
                eval(['mkdir ' foldername])

                %% select_intervals in LC-LC plane
               [ints,x,y,xmin,xmax,ymin,ymax] = select_interval(mz1,ind_mz1,xwindow,ywindow,threshold_peak);

                %% Evaluate and save results
                for ii = 1:length(ints)

                    ind_mz1 = find(mz1(x(ii),y(ii),:) ~= 0);
                    mz1_temp = mz1(xmin(ii):xmax(ii),ymin(ii):ymax(ii),ind_mz1);
                    axis_mz1_temp = axis_mz1(ind_mz1);
                    ind_mz2 = find(mz2(x(ii),y(ii),:) ~= 0);
                    mz2_temp = mz2(xmin(ii):xmax(ii),ymin(ii):ymax(ii),ind_mz2);
                    axis_mz2_temp = axis_mz2(ind_mz2);

                    [~,ind_target_1_temp] = min(abs(axis_mz1_temp-mz_target));
                    %  [~,ind_target_2_temp] = min(abs(axis_mz2_temp-mz_target));

                    ref_profile_lc_1 = sum(double(mz1_temp(:,:,ind_target_1_temp)),1);
                    ref_profile_lc_1 = ref_profile_lc_1./norm(ref_profile_lc_1,'fro');
                    ref_profile_lc_2 = sum(double(mz1_temp(:,:,ind_target_1_temp)),2);
                    ref_profile_lc_2 = ref_profile_lc_2./norm(ref_profile_lc_2,'fro');

                    %% MS1

                    test_profiles_lc2 = squeeze(sum(double(mz1_temp),2));
                    t_norm_lc2 = vecnorm(test_profiles_lc2,2,1);
                    test_profiles_lc2 = test_profiles_lc2./t_norm_lc2;

                    test_profiles_lc1 = squeeze(sum(double(mz1_temp),1));
                    t_norm_lc1 = vecnorm(test_profiles_lc1,2,1);
                    test_profiles_lc1 = test_profiles_lc1./t_norm_lc1;

                    coef_lc2_mz1 = ref_profile_lc_2'*test_profiles_lc2;
                    coef_lc1_mz1 = ref_profile_lc_1*test_profiles_lc1;

                    score_all_mz1 = coef_lc1_mz1.*coef_lc2_mz1;
                    masses_to_keep = find(score_all_mz1 > threshold_masses);

                    mz1_filtered = mz1_temp(:,:,masses_to_keep);
                    mz1_mass_spectrum = squeeze(sum(sum(double(mz1_filtered),1),2));
                    mz1_mass_spectrum = mz1_mass_spectrum./norm(mz1_mass_spectrum,'fro');
                    mz1_el_profil_lc1 = sum(sum(double(mz1_filtered),3),1);
                    mz1_el_profil_lc2 = sum(sum(double(mz1_filtered),3),2);
                    axis_mz1_filtered = axis_mz1_temp(masses_to_keep);

                    %% MS2

                    test_profiles_lc2 = squeeze(sum(double(mz2_temp),2));
                    t_norm_lc2 = vecnorm(test_profiles_lc2,2,1);
                    test_profiles_lc2 = test_profiles_lc2./t_norm_lc2;

                    test_profiles_lc1 = squeeze(sum(double(mz2_temp),1));
                    t_norm_lc1 = vecnorm(test_profiles_lc1,2,1);
                    test_profiles_lc1 = test_profiles_lc1./t_norm_lc1;

                    coef_lc2_mz2 = ref_profile_lc_2'*test_profiles_lc2;
                    coef_lc1_mz2 = ref_profile_lc_1*test_profiles_lc1;

                    score_all_mz2 = coef_lc1_mz2.*coef_lc2_mz2;
                    masses_to_keep = find(score_all_mz2 > threshold_masses);

                    mz2_filtered = mz2_temp(:,:,masses_to_keep);
                    mz2_mass_spectrum = squeeze(sum(sum(double(mz2_filtered),1),2));
                    mz2_mass_spectrum = mz2_mass_spectrum./norm(mz2_mass_spectrum,'fro');
                    mz2_el_profil_lc1 = sum(sum(double(mz2_filtered),3),1);
                    mz2_el_profil_lc2 = sum(sum(double(mz2_filtered),3),2);
                    axis_mz2_filtered = axis_mz2_temp(masses_to_keep);


                    %% saving the output as a structure file containing:
                    % - the reduced size, sparse tensors for MS1 and MS2 
                    %   defined by the selected intervals in LC-LC plane
                    % - the even more reduced size tensors after applying
                    %   mass filtering for MS1 and MS2
                    % - mz-spectra and mass axis and elution profiles (LC-1 / LC-2) after
                    %   mass filtering for MS1 and MS2

                    cd(foldername)

                    eval(['output' num2str(ii) '.mz1_temp = mz1_temp;'])
                    eval(['output' num2str(ii) '.mz2_temp = mz2_temp;'])
                    eval(['output' num2str(ii) '.axis_mz1_temp = axis_mz1_temp;'])
                    eval(['output' num2str(ii) '.axis_mz2_temp = axis_mz2_temp;'])

                    eval(['output' num2str(ii) '.mz1_filtered = mz1_filtered;'])
                    eval(['output' num2str(ii) '.mz2_filtered = mz2_filtered;'])
                    eval(['output' num2str(ii) '.axis_mz1_filtered = axis_mz1_filtered;'])
                    eval(['output' num2str(ii) '.axis_mz2_filtered = axis_mz2_filtered;'])

                    eval(['output' num2str(ii) '.mz1_mass_spectrum = mz1_mass_spectrum;'])
                    eval(['output' num2str(ii) '.mz2_mass_spectrum = mz2_mass_spectrum;'])
                    eval(['output' num2str(ii) '.mz1_el_profil_lc1 = mz1_el_profil_lc1;'])
                    eval(['output' num2str(ii) '.mz1_el_profil_lc2 = mz1_el_profil_lc2;'])
                    eval(['output' num2str(ii) '.mz2_el_profil_lc1 = mz2_el_profil_lc1;'])
                    eval(['output' num2str(ii) '.mz2_el_profil_lc2 = mz2_el_profil_lc2;'])

                    save(strcat('output', suspect_info_id{i}, num2str(ii)),strcat('output', num2str(ii)))

                    cd ..

                end


            end
        end
    
    end
end
