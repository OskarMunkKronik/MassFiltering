classdef DataAnalysis_class
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        MF_obj
        ID
        ind_data
        options_SIT
        options_MF
        BasePath
        SuspectList
    end

    methods
        function obj = init(obj,MF_obj)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.MF_obj = MF_obj;


            disp('please select "options_SIT.txt" file.')
            [filename,pathname] = uigetfile('*txt');
            filename = fullfile(pathname,filename);
            opt_table  = readtable(filename);

            obj.options_SIT.MaxIter = opt_table.MaxIter;
            obj.options_SIT.ConvCrit = opt_table.ConvCrit;
            obj.options_SIT.Constr = [opt_table.Constr_1,opt_table.Constr_2];
            obj.options_SIT.Init = opt_table.Init;
            obj.options_SIT.classify = opt_table.classify;
            if isfield(opt_table,'Fac')
                obj.options_SIT.Fac = opt_table.Fac;
            else
                obj.options_SIT.Fac = [];
            end
            clear opt_table

            disp('please select "options_MF.txt" file.')
            [filename,pathname] = uigetfile('*txt');
            filename = fullfile(pathname,filename);
            opt_table  = readtable(filename);

            obj.options_MF.threshold_peak = opt_table.threshold_peak;
            obj.options_MF.threshold_mass = opt_table.threshold_mass;
            obj.options_MF.threshold_sim_1 = opt_table.threshold_sim_1;
            obj.options_MF.threshold_sim_2 = opt_table.threshold_sim_2;
            obj.options_MF.xwindow = opt_table.xwindow;
            obj.options_MF.ywindow = opt_table.ywindow;

            clear opt_table

            disp('please select "suspect_list.xlsx" file.')
            [filename,pathname] = uigetfile('*xlsx');
            filename = fullfile(pathname,filename);
            suspect_list_raw = readtable(filename,'ReadVariableNames', true);

            %define column names for M/Z-Information
            column_name_mz = input('define column names for M/Z-Information: '); % uncomment this
            eval(strcat('suspect_info_masses = suspect_list_raw.',column_name_mz))

            %define column names for analyte names
            column_name_id = input('define column names for analyte names: '); % uncomment this
            eval(strcat('suspect_info_id = suspect_list_raw.',column_name_id))
            sanitizeString = @(str) regexprep(str,{,',',';','[',']','(',')',' ','+'},'');
            suspect_info_id = cellfun(sanitizeString,suspect_info_id,'UniformOutput',false);

            obj.SuspectList.suspect_info_id = suspect_info_id;
            obj.SuspectList.suspect_info_masses = suspect_info_masses;
            clear suspect_list_raw suspect_info_id column_name_id column_name_mz
            disp('please select the folder were to store the result files.')
            obj.BasePath = uigetdir();
            cd(obj.BasePath)
            disp('class has been successfully initialized.')
        end

        function obj = check_peaks(obj,screening)

            threshold_peak = obj.options_MF.threshold_peak;
            threshold_mass = obj.options_MF.threshold_mass;

            mz1 = obj.MF_obj.Reshaped_Data{1,1};
            axis_mz1 = obj.MF_obj.mz_temp{1,1};
            mz2 = obj.MF_obj.Reshaped_Data{1,2};
            axis_mz2 = obj.MF_obj.mz_temp{1,2};

            obj.options_SIT.peak_position    = zeros(size(obj.SuspectList.suspect_info_masses,1),2);
            obj.options_SIT.peak_position_rt = zeros(size(obj.SuspectList.suspect_info_masses,1),2);
             
            len_mz  = length(unique(regexprep(obj.SuspectList.suspect_info_id, '_peak_\d+$', '')));%size(obj.SuspectList.suspect_info_masses,1);
            for i = 1:len_mz

                mz_target = obj.SuspectList.suspect_info_masses(i);

                [dif_mz1,ind_mz1] = min(abs(axis_mz1-mz_target));

                if (dif_mz1/mz_target*10^6) >=  threshold_mass
%                     obj.options_SIT.peak_position(i,1) = 0;%[];
%                     obj.options_SIT.peak_position(i,2) = 0;% [];
                    display('suspect NOT found')
                else
                    %%
%                     [ints,x,y] = obj.peaks2(squeeze(double(mz1(:,:,ind_mz1))),'MinPeakHeight',threshold_peak);
  
                    %% Smooths data before peak detections
%                     clf
                    z  = smoothdata(squeeze(double(mz1(:,:,ind_mz1))),1,'gaussian',7);
                    [ints,x,y] = obj.peaks2(z,'MinPeakHeight',threshold_peak);

%%
%                     try
                        [~,inds] = sort(ints,'descend');
                        if length(inds)>20
                        inds = inds(1:20);
                        end 
                        ints = ints(inds);
                        x = x(inds);
                        y = y(inds);
                        if screening == 0
                            figure
                            mesh(squeeze(double(mz1(:,:,ind_mz1))))
                            hold on
                            ii = 1;
                            scatter3(y(ii),x(ii),ints(ii),'filled')
                            prompt_str = input("Enter 'y' for seein the next point and enter 'n' for moving on to next suspect: ","s")

                            if prompt_str == 'y'
                                while prompt_str == 'y' & ii <= length(y)
                                    ii = ii+1;
                                    hold on
                                    scatter3(y(ii),x(ii),ints(ii),'filled')
                                    prompt_str = input("Enter 'y' for seein the next point and enter 'n' for moving on to next suspect: ","s")
                                end
                            end


                            if ~isempty(x)
                                disp(i)
                                disp('suspect found')
                                obj.options_SIT.peak_position(i,1) = x(ii);
                                obj.options_SIT.peak_position(i,2) = y(ii);
                            else
                                disp(i)
                                disp('suspect NOT found')
                                obj.options_SIT.peak_position(i,1) = [];
                                obj.options_SIT.peak_position(i,2) = [];
                            end

                            x = [];
                            y = [];
                            ints = [];
                            close all

                        elseif screening == 1
                            
                            if ~isempty(x)
                                len = size(obj.SuspectList.suspect_info_masses,1);
                                len_x = length(x);
                                
                               if len_x>1
%                                    for ii = 1:len_x
                                       
%                                         if ii == 1
                                           cleanName = regexprep(obj.SuspectList.suspect_info_id{i}, '_peak_*\d+$', '');
%                                            Ind = find(ismember(obj.SuspectList.suspect_info_id,obj.SuspectList.suspect_info_id{i}),1,'first');
%                                            obj.SuspectList.suspect_info_id{i}      = [cleanName,'_peak_',num2str(1)];  
                                           
                                            newName = {[cleanName,'_peak_1']};
                                            Ind = find(ismember(obj.SuspectList.suspect_info_id,newName));
                                            obj.options_SIT.peak_position(i,1) = x(1);
                                           obj.options_SIT.peak_position(i,2) = y(1);
                                            if ~isempty(Ind)
                                                obj.SuspectList.suspect_info_id(Ind)  = newName;
                                            else

                                             obj.SuspectList.suspect_info_id(i)      = newName; %Fix this! 
                                            end 
%                                         else 
%                                         end 
%                                    end 
                                        for ii = 2:len_x
                                            
                                            newName = {[cleanName,'_peak_',num2str(ii)]};
                                            Ind = find(ismember(obj.SuspectList.suspect_info_id,newName));
                                            if ~isempty(Ind)
                                                obj.SuspectList.suspect_info_id(Ind)  = newName;
                                                obj.options_SIT.peak_position(Ind,:)       = [x(ii),y(ii)];
                                                obj.options_SIT.peak_position_rt(Ind,:)      = zeros(1,2);
                                                obj.SuspectList.suspect_info_masses(Ind)  = repelem(mz_target,1);
                                            else

                                             obj.SuspectList.suspect_info_id(len+ii-1)         = newName; %Fix this! 
                                             obj.options_SIT.peak_position(len+ii-1,:)         = [x(ii),y(ii)];
                                             obj.options_SIT.peak_position_rt(len+ii-1,:)      = zeros(1,2);
                                             obj.SuspectList.suspect_info_masses(len+ii-1)     = repelem(mz_target,1);
                                            end  
                                        end 
                                     
                                    
                                   
                               else 
                                disp(i)
                                disp('suspect found')
                                obj.options_SIT.peak_position(i,1) = x(1);
                                obj.options_SIT.peak_position(i,2) = y(1);
                                cleanName = regexprep(obj.SuspectList.suspect_info_id{i}, '_peak_*\d+$', '');
                                obj.SuspectList.suspect_info_id(i)      = {[cleanName,'_peak_1']};
                 
                             end 
                               
                                else
                                disp(i)
                                disp('suspect NOT found')
%                                 obj.options_SIT.peak_position(i,1) = 0;%[];
%                                 obj.options_SIT.peak_position(i,2) = 0;%[];
                            end

%                         end
                    end
                end
            end
            dropZeros                                       = sum(obj.options_SIT.peak_position,2)==0;
            obj.options_SIT.peak_position(dropZeros,:)      = [];
            obj.options_SIT.peak_position_rt(dropZeros,:)      = [];
            obj.SuspectList.suspect_info_masses(dropZeros)  = [];
            obj.SuspectList.suspect_info_id(dropZeros)      = [];

        end

        function obj = RunMassFiltering(obj,filtered,cutoff)

            threshold_peak = obj.options_MF.threshold_peak;

            if filtered == 0
                threshold_sim = obj.options_MF.threshold_sim_1;
            elseif filtered == 1
                threshold_sim = obj.options_MF.threshold_sim_2;
            end

            mz1 = obj.MF_obj.Reshaped_Data{1,1};
            axis_mz1 = obj.MF_obj.mz_temp{1,1};
            mz2 = obj.MF_obj.Reshaped_Data{1,2};
            axis_mz2 = obj.MF_obj.mz_temp{1,2};
            if ~isfield(obj.options_SIT,"peak_position")
                obj.options_SIT.peak_position = zeros(size(obj.SuspectList.suspect_info_masses,1),2);

                for i = 1:size(obj.SuspectList.suspect_info_masses,1)

                    mz_target = obj.SuspectList.suspect_info_masses(i);

                    [~,ind_mz1] = min(abs(axis_mz1-mz_target));
                    [~,ind_mz2] = min(abs(axis_mz2-mz_target));

                    [ints,x,y] = obj.peaks2(squeeze(double(mz1(:,:,ind_mz1))),'MinPeakHeight',threshold_peak);
                    [~,inds] = sort(ints,'descend');
                    ints = ints(inds);
                    x = x(inds);
                    y = y(inds);

                    obj.options_SIT.peak_position(i,1) = x(1);
                    obj.options_SIT.peak_position(i,2) = y(1);

                end
            end

            for i = 1:size(obj.SuspectList.suspect_info_masses,1)
               %Oskar edit cosme
                totalLength = size(obj.SuspectList.suspect_info_masses,1);  % Replace with the actual length
                  
                if mod(i, 10) == 0 || i == totalLength
                    fprintf('Processing %d / %d\n', i, totalLength);
                end
                % Your processing code here
               
                if ~isempty(obj.options_SIT.peak_position(i,1))
                    mz_target = obj.SuspectList.suspect_info_masses(i);
                    xmin = obj.options_SIT.peak_position(i,1)-obj.options_MF.xwindow;
                    ind_fix_x = find((xmin) < 1);
                    if ~isempty(ind_fix_x)
                        xmin(ind_fix_x) = 1;
                    end
                    xmax = obj.options_SIT.peak_position(i,1)+obj.options_MF.xwindow;
                    ind_fix_x = find((xmax) > size(mz1,1));
                    if ~isempty(ind_fix_x)
                        xmax(ind_fix_x) = size(mz1,1);
                    end

                    ymin = obj.options_SIT.peak_position(i,2)-obj.options_MF.ywindow;
                    ind_fix_y = find((ymin) < 1);
                    if ~isempty(ind_fix_y)
                        ymin(ind_fix_y) = 1;
                    end
                    ymax = obj.options_SIT.peak_position(i,2)+obj.options_MF.ywindow;
                    ind_fix_y = find((ymax) > size(mz1,1));
                    if ~isempty(ind_fix_y)
                        ymax(ind_fix_y) = size(mz1,1);
                    end

                    ind_mz1 = find(mz1(obj.options_SIT.peak_position(i,1),obj.options_SIT.peak_position(i,2),:) ~= 0);
                    mz1_temp = mz1(xmin:xmax,ymin:ymax,ind_mz1);
                    axis_mz1_temp = axis_mz1(ind_mz1);
                    ind_mz2 = find(mz2(obj.options_SIT.peak_position(i,1),obj.options_SIT.peak_position(i,2),:) ~= 0);
                    mz2_temp = mz2(xmin:xmax,ymin:ymax,ind_mz2);
                    axis_mz2_temp = axis_mz2(ind_mz2);

                    if ~isempty(cutoff)
                        mz1_temp = double(mz1_temp);
                        mz1_temp = mz1_temp(:,:,full(axis_mz1_temp) < mz_target+cutoff);
                        mz2_temp = double(mz2_temp);
                        mz2_temp = mz2_temp(:,:,full(axis_mz2_temp) < mz_target+cutoff);
                        axis_mz1_temp = axis_mz1_temp(axis_mz1_temp < mz_target+cutoff);
                        axis_mz2_temp = axis_mz2_temp(axis_mz2_temp < mz_target+cutoff);
                    end

                    [~,ind_target_1_temp] = min(abs(axis_mz1_temp-mz_target));

                    ref_profile_lc_1 = sum(double(mz1_temp(:,:,ind_target_1_temp)),1);
                    ref_profile_lc_1 = ref_profile_lc_1./norm(ref_profile_lc_1,'fro');
                    ref_profile_lc_2 = sum(double(mz1_temp(:,:,ind_target_1_temp)),2);
                    ref_profile_lc_2 = ref_profile_lc_2./norm(ref_profile_lc_2,'fro');

                    %% mz-1
                    test_profiles_lc2 = zeros(size(mz1_temp,[1,3]));
                    test_profiles_lc2(:,:) = squeeze(sum(double(mz1_temp),2));
                    t_norm_lc2 = vecnorm(test_profiles_lc2,2,1);
                    test_profiles_lc2 = test_profiles_lc2./t_norm_lc2;
                    
                    test_profiles_lc1 = zeros(size(mz1_temp,[2,3]));
                    test_profiles_lc1(:,:) = squeeze(sum(double(mz1_temp),1));
                    t_norm_lc1 = vecnorm(test_profiles_lc1,2,1);
                    test_profiles_lc1 = test_profiles_lc1./t_norm_lc1;

                    coef_lc2_mz1 = ref_profile_lc_2'*test_profiles_lc2;
                    coef_lc1_mz1 = ref_profile_lc_1*test_profiles_lc1;

                    score_all_mz1 = coef_lc1_mz1.*coef_lc2_mz1;
                    masses_to_keep = find(score_all_mz1 > threshold_sim);

                    [~,ind_lc1] = max(ref_profile_lc_1);
                    [~,ind_lc2] = max(ref_profile_lc_2);

                    mz1_filtered = mz1_temp(:,:,masses_to_keep);
                    mz1_mass_spectrum = double(squeeze(mz1_temp(ind_lc2,ind_lc1,masses_to_keep)));

                    mz1_mass_spectrum = mz1_mass_spectrum./max(mz1_mass_spectrum);
                    mz1_el_profil_lc1 = sum(sum(double(mz1_filtered),3),1);
                    mz1_el_profil_lc2 = sum(sum(double(mz1_filtered),3),2);
                    axis_mz1_filtered = axis_mz1_temp(masses_to_keep);

                    %% mz-2
                    test_profiles_lc2 = zeros(size(mz2_temp,[1,3]));
                    test_profiles_lc2(:,:) = squeeze(sum(double(mz2_temp),2));
                    t_norm_lc2 = vecnorm(test_profiles_lc2,2,1);
                    test_profiles_lc2 = test_profiles_lc2./t_norm_lc2;
                    
                    test_profiles_lc1 = zeros(size(mz2_temp,[2,3]));
                    test_profiles_lc1(:,:) = squeeze(sum(double(mz2_temp),1));
                    t_norm_lc1 = vecnorm(test_profiles_lc1,2,1);
                    test_profiles_lc1 = test_profiles_lc1./t_norm_lc1;

                    coef_lc2_mz2 = ref_profile_lc_2'*test_profiles_lc2;
                    coef_lc1_mz2 = ref_profile_lc_1*test_profiles_lc1;

                    score_all_mz2 = coef_lc1_mz2.*coef_lc2_mz2;

                    masses_to_keep = find(score_all_mz2 > threshold_sim);

                    mz2_filtered = mz2_temp(:,:,masses_to_keep);
                    mz2_mass_spectrum = double(squeeze(mz2_temp(ind_lc2,ind_lc1,masses_to_keep)));

                    mz2_mass_spectrum = mz2_mass_spectrum./max(mz2_mass_spectrum);
                    mz2_el_profil_lc1 = sum(sum(double(mz2_filtered),3),1);
                    mz2_el_profil_lc2 = sum(sum(double(mz2_filtered),3),2);
                    axis_mz2_filtered = axis_mz2_temp(masses_to_keep);

                    peak_position     = obj.options_SIT.peak_position(i,:);
                    peak_position_rt  = obj.options_SIT.peak_position_rt(i,:);
                    %% saving the output
                    cd(obj.BasePath)
                    mkdir(obj.SuspectList.suspect_info_id{i})
                    cd(obj.SuspectList.suspect_info_id{i})

                    eval('output_MF.mz1_temp = mz1_temp;')
                    eval('output_MF.mz2_temp = mz2_temp;')
                    eval('output_MF.axis_mz1_temp = axis_mz1_temp;')
                    eval('output_MF.axis_mz2_temp = axis_mz2_temp;')

                    eval('output_MF.mz1_filtered = mz1_filtered;')
                    eval('output_MF.mz2_filtered = mz2_filtered;')
                    eval('output_MF.axis_mz1_filtered = axis_mz1_filtered;')
                    eval('output_MF.axis_mz2_filtered = axis_mz2_filtered;')

                    eval('output_MF.mz1_mass_spectrum = mz1_mass_spectrum;')
                    eval('output_MF.mz2_mass_spectrum = mz2_mass_spectrum;')
                    eval('output_MF.mz1_el_profil_lc1 = mz1_el_profil_lc1;')
                    eval('output_MF.mz1_el_profil_lc2 = mz1_el_profil_lc2;')
                    eval('output_MF.mz2_el_profil_lc1 = mz2_el_profil_lc1;')
                    eval('output_MF.mz2_el_profil_lc2 = mz2_el_profil_lc2;')
                    eval('output_MF.mz_target = mz_target;')
                    eval('output_MF.ind_mz1 = ind_mz1;')
                    eval('output_MF.ind_mz2 = ind_mz2;')
                    
                    eval('output_MF.peak_position = peak_position;')
                    eval('output_MF.peak_position_rt = peak_position_rt;')
                    
                    if filtered == 0
                        save('output_MF','output_MF')
                    else
                        save('output_ignore_MF','output_MF')
                    end

                    cd ..
                end
            end
        end

        function obj = get_fac(obj,filtered)
            % I want to get the number of factors for each SIT model that
            % should be fitted
            cd(obj.BasePath);
            for i = 1:size(obj.SuspectList.suspect_info_masses,1)

                cd(obj.SuspectList.suspect_info_id{i})
                if filtered == 0
                    load('output_MF.mat')
                else
                    load('output_ignore_MF.mat')
                end
                if~isempty(output_MF.mz2_mass_spectrum)
                [Fac11,Fac12,Fac13,Fac21,Fac22,Fac23] = obj.get_rank_of_modes(output_MF,filtered);
                obj.options_SIT.Fac(1,i) = min(Fac11,Fac12);
                obj.options_SIT.Fac(2,i) = min(Fac21,Fac22);

                factors_of_modes = table(Fac11,Fac12,Fac13,Fac21,Fac22,Fac23);
                if filtered == 0
                    writetable(factors_of_modes,'factors_of_modes','FileType','spreadsheet')
                else
                    writetable(factors_of_modes,'factors_of_modes_MF_SIT','FileType','spreadsheet')
                end
                end
                cd ..

            end

        end

        function obj = run_mcr(obj,filtered,model_type)

            cd(obj.BasePath)

            modes = {'mz1';'mz2'};
            mz_target = [];
            options.MaxIter    = obj.options_SIT.MaxIter;
            options.ConvCrit   = obj.options_SIT.ConvCrit;
            options.Constr     = obj.options_SIT.Constr;
            options.Init       = obj.options_SIT.Init;
            options.classify   = obj.options_SIT.classify;
            options.InitLoads  = [];

            for i = 1:size(obj.SuspectList.suspect_info_masses,1)
                cd(obj.SuspectList.suspect_info_id{i})
                if filtered == 0
                    load('output_MF.mat');
                else
                    load('output_ignore_MF.mat');
                end
                if ~isempty(output_MF.mz2_mass_spectrum)
                for ii = 1:size(modes,1)
                    if filtered == 0
                        eval("data = output_MF.([modes{ii} '_temp']);")
                        eval("mz_axis = output_MF.(['axis_' modes{ii} '_temp']);")
                        %         eval("data = output1.([modes{ii} '_filtered']);")
                        %         eval("mz_axis = output1.(['axis_' modes{ii} '_filtered']);")
                        eval("mz_target = output_MF.mz_target;")
                    else
                        eval("data = output_MF.([modes{ii} '_filtered']);")
                        eval("mz_axis = output_MF.(['axis_' modes{ii} '_filtered']);")
                        %         eval("data = output1.([modes{ii} '_filtered']);")
                        %         eval("mz_axis = output1.(['axis_' modes{ii} '_filtered']);")
                        eval("mz_target = output_MF.mz_target;")
                    end


                    [~,ind_mz_target] = min(abs(mz_axis-mz_target));
                    data = double(data);
                    if ~(length(size(data)) < 3 | min(size(data)) < 2)
                        Fac = obj.options_SIT.Fac(ii,i);

                        if model_type == "MCR"

                            for iii = 1:5
                                model{iii} = SIT_MCR(data,Fac,options);
                                display(['compound: ' num2str(i) '/' num2str(size(obj.SuspectList.suspect_info_masses,1)) ' mode: ' num2str(ii) ' model: ' num2str(iii)])
                            end

                            cd(obj.BasePath)
                            cd(obj.SuspectList.suspect_info_id{i})
                        elseif model_type == "SIT"
                            for iii = 1:5
                                model{iii} = DataAnalysis_class.SIT(data,Fac,options);
                                display(['compound: ' num2str(i) '/' num2str(size(obj.SuspectList.suspect_info_masses,1)) ' mode: ' num2str(ii) ' model: ' num2str(iii)])
                            end
                            cd(obj.BasePath)
                            cd(obj.SuspectList.suspect_info_id{i})

                        end

                        [sp,Xhat] = obj.createOutputFromModel(model,ind_mz_target);

                        eval(['output_SIT.' modes{ii} '_mass_spectrum' '= sp;'])
                        eval(['output_SIT.' modes{ii} '_Xhat' '= Xhat;'])
                        eval(['output_SIT.' modes{ii} '_axis_temp = mz_axis;'])
                        eval(['output_SIT.mz_target = mz_target'])
                        clear model

                        if filtered == 0
                            save('output_SIT.mat','output_SIT');
                        else
                            save('output_MF_SIT.mat','output_SIT');
                        end
                    end 
                end
                       
                    end
                 cd ..
            end

        end

        function obj = get_results(obj)

            cd(obj.BasePath)
            for i = 1:size(obj.SuspectList.suspect_info_masses,1)

                cd(obj.SuspectList.suspect_info_id{i})

                if isfile('output_MF.mat')
                    mtype = 1;
                    load('output_MF.mat');
                    outputFile = 'spectra_MF.xlsx';
                    obj.make_result_tables(output_MF,outputFile,mtype)

                    outputFile = 'spectra_PAM.xlsx';
                    obj.PeakApexMethod(output_MF,outputFile)
                    clear outputFile output_MF
                end

                if isfile('output_SIT.mat')
                    mtype = 2;
                    load('output_SIT.mat');
                    outputFile = 'spectra_SIT.xlsx';
                    obj.make_result_tables(output_SIT,outputFile,mtype)
                    clear outputFile output_SIT

                end
                if isfile('output_MF_SIT.mat')
                    mtype = 2;
                    load('output_MF_SIT.mat');
                    outputFile = 'spectra_MF_SIT.xlsx';
                    obj.make_result_tables(output_SIT,outputFile,mtype)
                    clear outputFile output_SIT

                end

                cd ..

            end

        end

        function obj = compareWithRefSpec(obj)

            cd(obj.BasePath)

            for ii = 1:size(obj.SuspectList.suspect_info_id,1)
                cd(obj.SuspectList.suspect_info_id{ii})

                if isfile("ref_spectrum.xlsx")
                    ref_spectrum = readtable("ref_spectrum.xlsx",Sheet="RefSpec");
                    mz_dif = obj.options_MF.threshold_mass;

                    MF_spectrum = readtable("spectra_MF.xlsx",Sheet="Sheet1");
                    MF_SIT_spectrum = readtable("spectra_MF_SIT.xlsx",Sheet="Sheet1");
                    PAM_spectrum = readtable("spectra_PAM.xlsx",Sheet="Sheet1");
                    SIT_spectrum = readtable("spectra_SIT.xlsx",Sheet="Sheet1");

                    figure
                    %plot PAM spectrum and Ref_spectrum overlayed and mark the fragments in PAM
                    %spectrum that are matches in a different colour
                    subplot(2,2,1)
                    plot_spectra(PAM_spectrum,ref_spectrum,mz_dif)

                    subplot(2,2,2)
                    plot_spectra(SIT_spectrum,ref_spectrum,mz_dif)

                    subplot(2,2,3)
                    plot_spectra(MF_spectrum,ref_spectrum,mz_dif)

                    subplot(2,2,4)
                    plot_spectra(MF_SIT_spectrum,ref_spectrum,mz_dif)

                    saveas(gcf, 'compareWithref.svg', 'svg');
                    close all
                    cd ..
                else
                    cd ..
                end
            end
        end

        function obj = statistical_analysis(obj)

            cd(obj.BasePath)

            for ii = 1:size(obj.SuspectList.suspect_info_id,1)

                cd(obj.SuspectList.suspect_info_id{ii})
                if isfile("ref_spectrum.xlsx")
                    threshold_mass = obj.options_MF.threshold_mass;
                    [similarity,shared_fragments,all_shared_fragments] = obj.refspectra_analysis(threshold_mass);
                    varnames = {'PAM','SIT','MF','MF+SIT'};
                    similarity = table(similarity(:,1),similarity(:,2),similarity(:,3),similarity(:,4),'VariableNames',varnames);
                    shared_fragments = table(shared_fragments(:,1),shared_fragments(:,2),shared_fragments(:,3),shared_fragments(:,4),'VariableNames',varnames);
                    all_shared_fragments = cell2table(all_shared_fragments);
                    writetable(all_shared_fragments,'all_shared_fragments.xlsx');
                    writetable(similarity,'similarity.xlsx');
                    writetable(shared_fragments,'shared_fragments.xlsx');
                end
                cd ..

            end

        end

    end


    methods(Static)


        function model = SIT(X,Fac,options)

            % I/O
            %
            % I:
            % X(elutiontimes x samples x spectral chanels) = Low rank chromatographic data set
            % Fac = Number of latent variables to be fitted
            % options.MaxIter   :   maximum number of iterations (default = 1000)
            % options.ConvCrit  :   convergence criterion (default = 1e-09)
            % options.Constr    :   [Spectra, Elutionprofiles&Scores]
            %                       e.g [1 1] (default)
            %                       0 = unconstrained
            %                       1 = non-negativity
            % options.Init      :   Initialization
            %                   :   0 = random
            %                   :   1 = pure sample (default)
            %                   :   2 = custom (requires options.InitLoads)
            % options.InitLoads :   insert customized Loadings
            % options.classify  :   1 = Runs a DeepNeuralNetwork to classify
            %                           elutionprofiles as '0=weird, 1=peak, 2=cutoff, 3=baseline';
            % options.PeakshapeCorrection : If set to vals > 1, tri-linearity
            %                               constraint is relaxed.
            % options.compression : makes sense if size(X,i)*size(X,j) << size(X,k)
            % options.compression.do : 0 = no compression, 1 = compression
            % options.compression.basis: orthogonal basis for compression
            %
            %O:
            %model.spectra                      = (Fac x spectral chanels)
            %model.elutionprofiles              = (elutiontimes x Fac x samples)
            %model.scores                       = (Fac x samples)
            %model.detail.fit.X_sum_sq          = Total Sum of Squares;
            %model.detail.fit.res_sum_sq        = Residual Sum of Squares;
            %model.detail.fit.PercVar           = Explained Variance [%];
            %model.detail.fit.fitdif            = Final difference in Fit;
            %model.detail.lossfunc              = Loss Function values over all iterations;
            %model.detail.lossfunc_dif          = Fit difference over all iterations;
            %model.detail.residuals             = Residuals in dimensions of input data;
            %model.detail.iterations            = Number of Iterations
            %model.detail.time                  = Computation time;
            %model.detail.Profiles.PeakType     = PeakType
            %model.detail.Profiles.Niceness     = Niceness
            %model.detail.Profiles.PeakTypeHelp = '0 = weird, 1 = peak, 2 = cutoff, 3 = baseline';


            %% Unfold X
            X1 = X;
            mgc    = size(X1);

            X1 = reshape(X1,[mgc(1)*mgc(2) mgc(3)]);
            %'Mode 1 is elution time * sample
            %'Mode 2 is m/z fragments
            %% Arguments and option settings
            tic()
            ncomp           = Fac;
            if nargin < 3
                init        = 1;
                maxit       = 1000;
                constr      = [1 1];
                convcrit    = 1e-09;
                pccf        = 1;
            else
                if isfield(options,'Init')
                    init        = options.Init;
                else
                    init    = 0;
                end
                maxit       = options.MaxIter;
                if isfield(options,'Constr')
                    constr      = options.Constr;
                else
                    constr = [1 1];
                end
                convcrit    = options.ConvCrit;
                if isfield(options,'PeakshapeCorrection')
                    pccf        = options.PeakshapeCorrection;
                else
                    pccf = 1;
                end
                if isfield(options,'compression')
                    compression = options.compression.do;
                    basis_f     = options.compression.basis;
                else
                    compression = 0;
                end

                if ~isfield(options,'classify');
                    options.classify = 1;
                end
            end

            %% Initialization of S or C
            % pure sample / spectra

            if init == 0
                CB = randn(size(X1,1),ncomp);
            elseif init == 1
                [CB0,ind]=pure(X1,ncomp,10);
                CB = CB0';
                [S0,ind]=pure(X1',ncomp,10);
                S = S0';
            elseif init == 2
                S = options.InitLoads;
                StS = S'*S;
                StXt = S'*X1';
                CB = fcnnls([],[],StS,StXt)';
            end

            SST = sum(X1.^2,'all');
            SSE = zeros(1,2);
            SSE(1) = 2*SST;
            SSE(2) = SST;
            fitdif = SSE(1)-SSE(2);
            fit = [];
            iter = [];
            tictic = 0;

            %% calculating mcr solution
            % X = CS'
            % S has Norm 1
            % S = S/||S||
            % X*S*pinv(S'S) = C --> pinv(S'S)*S'X' = C;
            % pinv(C'C)*C'*X = S --> pinv(C'C)*C'*X = S;
            it = 0;
            while fitdif > convcrit & it < maxit
                it = it+1;
                CtC = CB'*CB;
                CtX = CB'*X1;
                if      constr(1) == 0
                    S = CtX'*pinv(CtC);
                    S(S<0) = 0;
                elseif  constr(1) == 1
                    if compression == 0
                        try
                            S = fcnnls(CB, X1, CtC, CtX)';
                            for k =1:ncomp
                                if sum(S(:,k)) == 0
                                    S(:,k) = rand(1,size(S,1));
                                end
                            end
                        catch
                            CB = randn(size(CB));
                            CtC = CB'*CB;
                            CtX = CB'*X1;
                            S = fcnnls(CB, X1, CtC, CtX)';
                        end
                    elseif compression == 1
                        E = CB;
                        f = X1;
                        for jj = 1:Fac
                            S(jj,:) = LSI(E,f(:,jj),basis_f,0);
                        end
                    end
                    for k = 1:ncomp
                        S(:,k) = S(:,k)*1/norm(S(:,k),'fro');
                    end
                end
                StS = S'*S;
                StXt = S'*X1';
                if      constr(2) == 0
                    CB = StXt'*pinv(StS);
                    CB(CB<0) = 0;
                elseif  constr(2) == 1
                    try
                        CB = fcnnls([],[],StS,StXt)';
                    catch
                        S = randn(size(S));
                        for k = 1:ncomp
                            S(:,k) = S(:,k)*1/norm(S(:,k),'fro');
                        end
                        StS = S'*S;
                        StXt = S'*X1';
                        CB = fcnnls([],[],StS,StXt)';
                    end
                end


                %% Shift invariant trilinearity on CB
                copt1 = [];
                Bft = [];
                Cft = [];
                Bneu = {};
                Cneu = [];

                for j = 1:ncomp;

                    [y,ya,p]= shiftmap(reshape(CB(:,j),mgc(1),mgc(2)));%shiftmap(copt1');
                    %making trilinearity for shifted data

                    [U,T,V] = svd(y,'econ');

                    %         if fitdif(end) < 100*convcrit | tictic > 0
                    %             tictic = 1;
                    %             ttemp = diag(T);
                    %             indstemp = find(ttemp./sum(ttemp) > 0.01);
                    %             if 1 > max(indstemp);
                    %                 pccf = 1;
                    %             else
                    %                 pccf = max(indstemp);
                    %             end
                    %         else
                    %             pccf = 1;
                    %         end
                    %         try
                    %             dif1(j,it) = pccf;
                    %         end
                    %trying flexibilisation of SIT
                    %                     ttemp = diag(T);
                    %                     indstemp = find(ttemp./sum(ttemp) > 0.01);
                    %                     if 1 > max(indstemp);
                    %                         pccf = 1;
                    %                     else
                    %                     pccf = max(indstemp);
                    %                     end
                    %         dif1(j,it) = pccf;
                    %         pccf = 1;
                    %storing results as elutionprofiles and relative concentrations
                    Bft = U(:,1:pccf)*T(1:pccf,1:pccf);
                    Cft = V(:,1:pccf)';
                    %shifting back
                    Bneu{1,j} = shiftmap(Bft*Cft,ya);
                    Bneu{1,j} = Bneu{1,j}(1:mgc(1),:);
                end
                %reshaping to get CB back
                BB = zeros(mgc(1),mgc(2), ncomp);
                for i = 1:(ncomp)
                    BB(:,:,i) = Bneu{1,i};
                end
                %           zeroforcing !!!! Spielwiese !!!!
                CBneu = reshape(BB,mgc(1)*mgc(2),ncomp);
                CB = CBneu;


                %evaluating loss function
                SSE(2) = SSE(1);
                %SSE(1) = sum((X1-CB*S').^2,'all');
                SSE(1) = SST+sum(sum((CB'*CB) .*(S'*S)))-2*sum((X1*S) .* CB,'all');
                fit(it) = (1-SSE(1)/SST)*100;
                fitdif(it) = abs((1-SSE(2)/SST)*100-(1-SSE(1)/SST)*100);

            end
            %Outputs
            % StS = S'*S;
            % StXt = S'*X1';
            % CB = Alg.fcnnls([],[],StS,StXt)';
            CB = reshape(CB,mgc(1),mgc(2),ncomp);
            model.spectra   = S;
            for i = 1:size(CB ,3)
                for ii = 1:size(CB ,2)
                    Cneu(ii,i) = norm(squeeze(CB(:,ii,i)),'fro');
                    CB(:,ii,i) = CB(:,ii,i)./Cneu(ii,i);
                end
            end
            model.elutionprofiles  = CB;
            model.elutionprofiles  = permute(model.elutionprofiles,[1 3 2]);
            % for i = 1:size(BB,3)
            %     for ii = 1:size(BB,2)
            %         Cneu(ii,i) = norm(squeeze(BB(:,ii,i)),'fro');
            %     end
            % end
            model.scores   = Cneu;

            model.detail.fit.X_sum_sq        = SST;
            model.detail.fit.res_sum_sq      = SSE(1);
            model.detail.fit.PercVar         = fit(end);
            model.detail.fit.fitdif          = fitdif(end);
            % model.detail.dif1                = dif1;
            % model.detail.lossfunc          = fit;
            % model.detail.lossfunc_dif      = fitdif;
            %         % !!! Spielwiese
            % model.detail.dif1              = dif1;
            % model.detail.dif2              = dif2;
            % Xhat                = makeXfromABC(model.spectra,model.elutionprofiles,model.scores);
            % model.detail.residuals     = permute(X,[2 1 3])-Xhat;

            model.detail.iterations    = it;
            model.detail.time          = toc()
            if options.classify == 1
                load('DeepNet2018b_1.mat');
                [~,PeakType,Niceness] = assessprof(model.elutionprofiles ,FinalDeepNet);
                model.detail.Profiles.PeakType = PeakType;
                model.detail.Profiles.Niceness = Niceness;
                model.detail.Profiles.PeakTypeHelp = '0 = weird, 1 = peak, 2 = cutoff, 3 = baseline';
            end





            %%


            function [y,ya,p] = shiftmap(x,p)
                %
                %  For calculating |fft| and phase map:
                %  INPUTS:
                %     x = MxN operates on the columns
                %  OPTIONAL INPUT:
                %     p = nonnegative scalar length of padding for fft
                %         if not included p = 2.^nextpow2(M);
                %  OUTPUTS:
                %     y = |fft(x)|
                %    ya = phase map of FFT: fft(x)./|fft(x)|
                %     p = p = 2.^nextpow2(n); [used when (p) not input].
                %
                %I/O: [y,ya,p] = shiftmap(x,p);
                %
                %  For calculating shifted profiles from |fft| and phase map:
                %  INPUTS:
                %     y = |fft(xhat)|, e.g., Xhat from a 1 PC PCA model of |FFT(x)|
                %    ya = phase map of FFT: fft(x)./|fft(x)|
                %  OUTPUTS:
                %     x = real(ifft(y.*ya))
                %
                %I/O: [x] = shiftmap(y,ya);

                yr    = 1e-6;             %regularization for angle map

                m     = size(x);

                if nargin<2||isscalar(p)  %calculate |FFT| and angle map
                    if nargin<2
                        p   = 2.^nextpow2(m(1));
                    end

                    % FFT
                    %             if isdataset(x)
                    %                 z   = fft(x.data,p);
                    %             else
                    %     z   = fft(x,p);       %assume x is class double
                    z = fft(x);
                    %             end
                    % |FFT| and angle map
                    y     = abs(z);
                    a     = y;
                    a(y<yr)   = yr;
                    ya    = z./a;
                else %p = ya (phase map)

                    % IFFT
                    %             if isdataset(x)
                    %                 y   = real(ifft(x.data.*p));
                    %             else
                    y   = real(ifft(x.*p));   %assume x is class double
                    %            end
                    ya    = [];
                    p     = [];


                end
            end

            function [sp,imp]=pure(d,nr,f)
                % [sp,imp]=pure(d,nr,f)
                % sp purest row/column profiles
                % imp indexes of purest variables
                % d data matrix; nr (rank) number of pure components to search
                % if d(nspectra,nwave) imp gives purest nwave => sp are conc. profiles (nr,nspectra)
                % if d(nwave,nspectra) imp gives purest nspectra => sp are spectra profiles (nr,nwave)
                % f percent of noise allowed respect maximum of the average spectrum given in % (i.e. 1% or 0.1%))

                [nrow,ncol]=size(d);


                % calculation of the purity spectrum

                f=f/100;
                s=std(d);
                m=mean(d);
                ll=s.*s+m.*m;
                f=max(m)*f;
                p=s./(m+f);

                [mp,imp(1)]=max(p);

                % calculation of the correlation matrix
                % l=sqrt(m.*m+(s+f).*(s+f));

                l=sqrt((s.*s+(m+f).*(m+f)));
                % dl=d./(l'*ones(1,ncol));
                for j=1:ncol,
                    dl(:,j)=d(:,j)./l(j);
                end
                c=(dl'*dl)./nrow;

                % calculation of the weights
                % first weight

                w(1,:)=ll./(l.*l);
                p(1,:)=w(1,:).*p(1,:);
                s(1,:)=w(1,:).*s(1,:);
                % figure(1)
                % subplot(3,1,1),plot(m)
                % title('unweigthed mean, std and first pure spectrum')
                % subplot(3,1,2),plot(s)
                % subplot(3,1,3),plot(p)

                % pause

                % next weights


                for i=2:nr,

                    for j=1:ncol,
                        [dm]=wmat(c,imp,i,j);
                        w(i,j)=det(dm);
                        p(i,j)=p(1,j).*w(i,j);
                        s(i,j)=s(1,j).*w(i,j);
                    end

                    % next purest and standard deviation spectrum

                    % plot(p(i,:))
                    % plot(s(i,:))
                    % title('sd and purest spectrum')
                    % figure(i)
                    % subplot(2,1,1),plot(p(i,:))
                    % title('next pure spectrum and std dev. spectrum')
                    % subplot(2,1,2),plot(s(i,:))
                    % pause


                    [mp(i),imp(i)]=max(p(i,:));
                end


                for i=1:nr,
                    impi=imp(i);
                    sp(1:nrow,i)=d(1:nrow,impi);
                end

                % figure(nr+1)
                sp=normv2(sp');
                % plot(sp')

            end

            function [dm]=wmat(c,imp,irank,jvar)
                dm(1,1)=c(jvar,jvar);
                for k=2:irank,
                    kvar=imp(k-1);
                    dm(1,k)=c(jvar,kvar);
                    dm(k,1)=c(kvar,jvar);
                    for kk=2:irank,
                        kkvar=imp(kk-1);
                        dm(k,kk)=c(kvar,kkvar);
                    end
                end
            end

            function [sn]=normv2(s)
                % normalitzacio s=s/sqrt(sum(si)2))
                [m,n]=size(s);
                for i=1:m,
                    sr=sqrt(sum(s(i,:).*s(i,:)));
                    sn(i,:)=s(i,:)./sr;
                end
            end


            % ****************************** Subroutine****************************

            function [X] = makeXfromABC(A,B,C)

                na = size(A,1);
                nb = size(B,1);
                nc = size(C,1);
                ncomp = size(C,2);

                X = zeros(na,nb,nc);
                for i = 1:nc
                    Di = zeros(ncomp,ncomp);
                    for ii = 1:ncomp
                        Di(ii,ii) = C(i,ii);
                    end
                    X(:,:,i) = A*Di*B(:,:,i)';
                end

                X = permute(X,[3,2,1]);
            end

            function x=LSI(E,f,G,h);

                % x=LSI(E,f,G,h);
                % Solves problem LSI min(norm(Ex-f))) subject to Gx>h
                %
                % From Lawson & Hanson 74, p. 167
                %
                % Copyright 1998
                % Rasmus Bro
                % rasmus@optimax.dk/rb@kvl.dk

                [I,J]=size(E);
                r=rank(E);
                [u,s,v]=svd(E,0);

                K1=v(:,1:r);
                Q1=u(:,1:r);
                R=s(1:r,1:r);
                f1=Q1'*f;
                GKR=G*K1*pinv(R);
                z=ldp(GKR,h-GKR*f1);

                %z=Ry-f1;

                y=pinv(R)*(z+f1);
                x=K1*y;
                function b=ldp(G,h);

                    % find min||b|| subject to Gb=>h
                    %
                    % From Lawson & Hanson 74, p. 165
                    %
                    % Copyright 1998
                    % Rasmus Bro
                    % rasmus@optimax.dk/rb@kvl.dk


                    [I,J]=size(G);
                    E=[G';h'];
                    f=[zeros(J,1);1];
                    uhat=fastnnls(E'*E,E'*f);
                    r=E*uhat-f;

                    if norm(r)==0
                        disp(' No solution to LDP problem')
                    elseif r(J+1)~=0
                        b=-r(1:J)/r(J+1);
                    else
                        b=NaN;
                    end
                end
            end

            function [K, Pset] = fcnnls(C, A, CtC, CtA)
                % NNLS using normal equations and the fast combinatorial strategy
                %
                % I/O: [K, Pset] = fcnnls(C, A);
                % K = fcnnls(C, A);
                %
                % C is the nObs x lVar coefficient matrix
                % A is the nObs x pRHS matrix of observations
                % K is the lVar x pRHS solution matrix
                % Pset is the lVar x pRHS passive set logical array
                %
                % Pset: set of passive sets, one for each column
                % Fset: set of column indices for solutions that have not yet converged
                % Hset: set of column indices for currently infeasible solutions
                % Jset: working set of column indices for currently optimal solutions
                %
                % Implementation is based on [1] with bugfixes, direct passing of sufficient stats,
                % and preserving the active set over function calls.
                %
                % [1] Van Benthem, M. H., & Keenan, M. R. (2004). Fast algorithm for the
                %   solution of large‐scale non‐negativity‐constrained least squares problems.
                %   Journal of Chemometrics: A Journal of the Chemometrics Society, 18(10), 441-450.


                % Check the input arguments for consistency and initialize
                if nargin == 2
                    error(nargchk(2,2,nargin))
                    [nObs, lVar] = size(C);

                    if size(A,1)~= nObs, error('C and A have imcompatible sizes'), end
                    if size(C,1) == size(C,2)
                        %         warning('A square matrix "C" was input, ensure this is on purpose.')
                    end
                    pRHS = size(A,2);
                    % Precompute parts of pseudoinverse
                    CtC = C'*C; CtA = C'*A;
                else
                    [lVar,pRHS] = size(CtA);

                end

                if nargin == 2 && size(C,1) == size(C,2)
                    warning('fcnnls - The coefficient matrix (C) was square - is this true or are you passing C''C?')
                end

                W = zeros(lVar, pRHS);
                iter = 0;
                maxiter = 6*lVar;


                % Obtain the initial feasible solution and corresponding passive set
                K = cssls(CtC, CtA);
                Pset=K>0;
                K(~Pset) = 0;
                D=K;
                Fset = find(~all(Pset));

                % Active set algorithm for NNLS main loop
                iter_outer = 1;
                while ~isempty(Fset) && iter_outer < maxiter
                    iter_outer = iter_outer + 1;
                    % Solve for the passive variables (uses subroutine below)
                    K(:,Fset) = cssls(CtC, CtA(:,Fset), Pset(:,Fset));
                    % Find any infeasible solutions
                    %     Hset = Fset(find(any(K(:,Fset) < 0)));
                    Hset = Fset((any(K(:,Fset) < 0)));

                    % Make infeasible solutions feasible (standard NNLS inner loop)
                    if ~isempty(Hset)
                        nHset = length(Hset);
                        alpha = zeros(lVar, nHset);

                        while ~isempty(Hset) && (iter < maxiter)
                            iter = iter + 1;
                            alpha(:,1:nHset) = Inf;
                            % Find indices of negative variables in passive set
                            [i, j] = find(Pset(:,Hset) & (K(:,Hset) < 0));
                            hIdx = sub2ind([lVar nHset], i, j);
                            %             if length(i) ~= length(j)
                            %                 keyboard
                            %             end
                            %             negIdx = sub2ind(size(K), i, Hset(j)'); % org
                            negIdx = sub2ind(size(K), i, reshape(Hset(j),size(i)));  % jlh mod
                            alpha(hIdx) = D(negIdx)./(D(negIdx) - K(negIdx));
                            [alphaMin,minIdx] = min(alpha(:,1:nHset));
                            alpha(:,1:nHset) = repmat(alphaMin, lVar, 1);
                            D(:,Hset) = D(:,Hset)-alpha(:,1:nHset).*(D(:,Hset)-K(:,Hset));
                            idx2zero = sub2ind(size(D), minIdx, Hset);
                            D(idx2zero) = 0;
                            Pset(idx2zero) = 0;
                            K(:, Hset) = cssls(CtC, CtA(:,Hset), Pset(:,Hset));
                            Hset = find(any(K < 0));
                            nHset = length(Hset);
                        end
                    end
                    % Make sure the solution has converged
                    %if iter == maxiter, warning('Maximum number iterations exceeded'), end
                    % Check solutions for optimality
                    W(:,Fset) = CtA(:,Fset)-CtC*K(:,Fset);
                    Jset = find(all(~Pset(:,Fset).*W(:,Fset) <= 0));
                    Fset = setdiff(Fset, Fset(Jset));
                    % For non-optimal solutions, add the appropriate variable to Pset
                    if ~isempty(Fset)
                        [mx, mxidx] = max(~Pset(:,Fset).*W(:,Fset));
                        Pset(sub2ind([lVar pRHS], mxidx, Fset)) = 1;
                        D(:,Fset) = K(:,Fset);
                    end
                end

                % ****************************** Subroutine****************************
                function [K] = cssls(CtC, CtA, Pset)
                    % Solve the set of equations CtA = CtC*K for the variables in set Pset
                    % using the fast combinatorial approach
                    K = zeros(size(CtA));
                    if (nargin == 2) || isempty(Pset) || all(Pset(:))
                        K = CtC\CtA; % Not advisable if matrix is close to singular or badly scaled
                        %     K = pinv(CtC)*CtA;
                    else
                        [lVar, pRHS] = size(Pset);
                        codedPset = 2.^(lVar-1:-1:0)*Pset;
                        [sortedPset, sortedEset] = sort(codedPset);
                        breaks = diff(sortedPset);
                        breakIdx = [0 find(breaks) pRHS];
                        for k = 1:length(breakIdx)-1
                            cols2solve = sortedEset(breakIdx(k)+1:breakIdx(k+1));
                            vars = Pset(:,sortedEset(breakIdx(k)+1));
                            K(vars,cols2solve) = CtC(vars,vars)\CtA(vars,cols2solve);
                            %         K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
                        end
                    end
                end
            end
        end

        function model = SIT_MCR(X,Fac,options)
            % I/O
            %
            % I:
            % X(elutiontimes x samples x spectral chanels) = Low rank chromatographic data set
            % Fac = Number of latent variables to be fitted
            % options.MaxIter   :   maximum number of iterations (default = 1000)
            % options.ConvCrit  :   convergence criterion (default = 1e-09)
            % options.Constr    :   [Spectra, Elutionprofiles&Scores]
            %                       e.g [1 1] (default)
            %                       0 = unconstrained
            %                       1 = non-negativity
            % options.Init      :   Initialization
            %                   :   0 = random
            %                   :   1 = pure sample (default)
            %                   :   2 = custom (requires options.InitLoads)
            % options.InitLoads :   insert customized Loadings
            % options.classify  :   1 = Runs a DeepNeuralNetwork to classify
            %                           elutionprofiles as '0=weird, 1=peak, 2=cutoff, 3=baseline';
            % options.PeakshapeCorrection : If set to vals > 1, tri-linearity
            %                               constraint is relaxed.
            % options.compression : makes sense if size(X,i)*size(X,j) << size(X,k)
            % options.compression.do : 0 = no compression, 1 = compression
            % options.compression.basis: orthogonal basis for compression
            %
            %O:
            %model.spectra                      = (Fac x spectral chanels)
            %model.elutionprofiles              = (elutiontimes x Fac x samples)
            %model.scores                       = (Fac x samples)
            %model.detail.fit.X_sum_sq          = Total Sum of Squares;
            %model.detail.fit.res_sum_sq        = Residual Sum of Squares;
            %model.detail.fit.PercVar           = Explained Variance [%];
            %model.detail.fit.fitdif            = Final difference in Fit;
            %model.detail.lossfunc              = Loss Function values over all iterations;
            %model.detail.lossfunc_dif          = Fit difference over all iterations;
            %model.detail.residuals             = Residuals in dimensions of input data;
            %model.detail.iterations            = Number of Iterations
            %model.detail.time                  = Computation time;
            %model.detail.Profiles.PeakType     = PeakType
            %model.detail.Profiles.Niceness     = Niceness
            %model.detail.Profiles.PeakTypeHelp = '0 = weird, 1 = peak, 2 = cutoff, 3 = baseline';


            %% Unfold X
            X1 = X;
            mgc    = size(X1);

            X1 = reshape(X1,[mgc(1)*mgc(2) mgc(3)]);
            %'Mode 1 is elution time * sample
            %'Mode 2 is m/z fragments
            %% Arguments and option settings
            tic()
            ncomp           = Fac;
            if nargin < 3
                init        = 1;
                maxit       = 1000;
                constr      = [1 1];
                convcrit    = 1e-09;
                pccf        = 1;
            else
                if isfield(options,'Init')
                    init        = options.Init;
                else
                    init    = 0;
                end
                maxit       = options.MaxIter;
                if isfield(options,'Constr')
                    constr      = options.Constr;
                else
                    constr = [1 1];
                end
                convcrit    = options.ConvCrit;
                if isfield(options,'PeakshapeCorrection')
                    pccf        = options.PeakshapeCorrection;
                else
                    pccf = 1;
                end
                if isfield(options,'compression')
                    compression = options.compression.do;
                    basis_f     = options.compression.basis;
                else
                    compression = 0;
                end

                if ~isfield(options,'classify');
                    options.classify = 1;
                end
            end

            %% Initialization of S or C
            % pure sample / spectra

            if init == 0
                CB = randn(size(X1,1),ncomp);
            elseif init == 1
                [CB0,ind]=pure(X1,ncomp,10);
                CB = CB0';
                [S0,ind]=pure(X1',ncomp,10);
                S = S0';
            elseif init == 2
                S = options.InitLoads;
                StS = S'*S;
                StXt = S'*X1';
                CB = fcnnls([],[],StS,StXt)';
            end

            SST = sum(X1.^2,'all');
            SSE = zeros(1,2);
            SSE(1) = 2*SST;
            SSE(2) = SST;
            fitdif = SSE(1)-SSE(2);
            fit = [];
            iter = [];
            tictic = 0;

            %% calculating mcr solution
            % X = CS'
            % S has Norm 1
            % S = S/||S||
            % X*S*pinv(S'S) = C --> pinv(S'S)*S'X' = C;
            % pinv(C'C)*C'*X = S --> pinv(C'C)*C'*X = S;
            it = 0;
            while fitdif > convcrit & it < maxit
                it = it+1;
                CtC = CB'*CB;
                CtX = CB'*X1;
                if      constr(1) == 0
                    S = CtX'*pinv(CtC);
                    S(S<0) = 0;
                elseif  constr(1) == 1
                    if compression == 0
                        try
                            S = fcnnls(CB, X1, CtC, CtX)';
                            for k =1:ncomp
                                if sum(S(:,k)) == 0
                                    S(:,k) = rand(1,size(S,1));
                                end
                            end
                        catch
                            CB = randn(size(CB));
                            CtC = CB'*CB;
                            CtX = CB'*X1;
                            S = fcnnls(CB, X1, CtC, CtX)';
                        end
                    elseif compression == 1
                        E = CB;
                        f = X1;
                        for jj = 1:Fac
                            S(jj,:) = LSI(E,f(:,jj),basis_f,0);
                        end
                    end
                    for k = 1:ncomp
                        S(:,k) = S(:,k)*1/norm(S(:,k),'fro');
                    end
                end
                StS = S'*S;
                StXt = S'*X1';
                if      constr(2) == 0
                    CB = StXt'*pinv(StS);
                    CB(CB<0) = 0;
                elseif  constr(2) == 1
                    try
                        CB = fcnnls([],[],StS,StXt)';
                    catch
                        S = randn(size(S));
                        for k = 1:ncomp
                            S(:,k) = S(:,k)*1/norm(S(:,k),'fro');
                        end
                        StS = S'*S;
                        StXt = S'*X1';
                        CB = fcnnls([],[],StS,StXt)';
                    end
                end


                %% Shift invariant trilinearity on CB
                copt1 = [];
                Bft = [];
                Cft = [];
                Bneu = {};
                Cneu = [];

                for j = 1:ncomp;

                    y= reshape(CB(:,j),mgc(1),mgc(2));%shiftmap(copt1');
                    %making trilinearity for shifted data

                    [U,T,V] = svd(y,'econ');

                    %
                    %storing results as elutionprofiles and relative concentrations
                    Bft = U(:,1:pccf)*T(1:pccf,1:pccf);
                    Cft = V(:,1:pccf)';
                    %shifting back
                    Bneu{1,j} = Bft*Cft;
                    Bneu{1,j} = Bneu{1,j}(1:mgc(1),:);
                end
                %reshaping to get CB back
                BB = zeros(mgc(1),mgc(2), ncomp);
                for i = 1:(ncomp)
                    BB(:,:,i) = Bneu{1,i};
                end
                %           zeroforcing !!!! Spielwiese !!!!
                CBneu = reshape(BB,mgc(1)*mgc(2),ncomp);
                CB = CBneu;


                %evaluating loss function
                SSE(2) = SSE(1);
                %SSE(1) = sum((X1-CB*S').^2,'all');
                SSE(1) = SST+sum(sum((CB'*CB) .*(S'*S)))-2*sum((X1*S) .* CB,'all');
                fit(it) = (1-SSE(1)/SST)*100;
                fitdif(it) = abs((1-SSE(2)/SST)*100-(1-SSE(1)/SST)*100);

            end
            %Outputs
            % StS = S'*S;
            % StXt = S'*X1';
            % CB = Alg.fcnnls([],[],StS,StXt)';
            CB = reshape(CB,mgc(1),mgc(2),ncomp);
            model.spectra   = S;
            for i = 1:size(CB ,3)
                for ii = 1:size(CB ,2)
                    Cneu(ii,i) = norm(squeeze(CB(:,ii,i)),'fro');
                    CB(:,ii,i) = CB(:,ii,i)./Cneu(ii,i);
                end
            end
            model.elutionprofiles  = CB;
            model.elutionprofiles  = permute(model.elutionprofiles,[1 3 2]);
            % for i = 1:size(BB,3)
            %     for ii = 1:size(BB,2)
            %         Cneu(ii,i) = norm(squeeze(BB(:,ii,i)),'fro');
            %     end
            % end
            model.scores   = Cneu;

            model.detail.fit.X_sum_sq        = SST;
            model.detail.fit.res_sum_sq      = SSE(1);
            model.detail.fit.PercVar         = fit(end);
            model.detail.fit.fitdif          = fitdif(end);
            % model.detail.dif1                = dif1;
            % model.detail.lossfunc          = fit;
            % model.detail.lossfunc_dif      = fitdif;
            %         % !!! Spielwiese
            % model.detail.dif1              = dif1;
            % model.detail.dif2              = dif2;
            % Xhat                = makeXfromABC(model.spectra,model.elutionprofiles,model.scores);
            % model.detail.residuals     = permute(X,[2 1 3])-Xhat;

            model.detail.iterations    = it;
            model.detail.time          = toc()
            if options.classify == 1
                load('DeepNet2018b_1.mat');
                [~,PeakType,Niceness] = assessprof(model.elutionprofiles ,FinalDeepNet);
                model.detail.Profiles.PeakType = PeakType;
                model.detail.Profiles.Niceness = Niceness;
                model.detail.Profiles.PeakTypeHelp = '0 = weird, 1 = peak, 2 = cutoff, 3 = baseline';
            end





            %%


            function [y,ya,p] = shiftmap(x,p)
                %
                %  For calculating |fft| and phase map:
                %  INPUTS:
                %     x = MxN operates on the columns
                %  OPTIONAL INPUT:
                %     p = nonnegative scalar length of padding for fft
                %         if not included p = 2.^nextpow2(M);
                %  OUTPUTS:
                %     y = |fft(x)|
                %    ya = phase map of FFT: fft(x)./|fft(x)|
                %     p = p = 2.^nextpow2(n); [used when (p) not input].
                %
                %I/O: [y,ya,p] = shiftmap(x,p);
                %
                %  For calculating shifted profiles from |fft| and phase map:
                %  INPUTS:
                %     y = |fft(xhat)|, e.g., Xhat from a 1 PC PCA model of |FFT(x)|
                %    ya = phase map of FFT: fft(x)./|fft(x)|
                %  OUTPUTS:
                %     x = real(ifft(y.*ya))
                %
                %I/O: [x] = shiftmap(y,ya);

                yr    = 1e-6;             %regularization for angle map

                m     = size(x);

                if nargin<2||isscalar(p)  %calculate |FFT| and angle map
                    if nargin<2
                        p   = 2.^nextpow2(m(1));
                    end

                    % FFT
                    %             if isdataset(x)
                    %                 z   = fft(x.data,p);
                    %             else
                    %     z   = fft(x,p);       %assume x is class double
                    z = fft(x);
                    %             end
                    % |FFT| and angle map
                    y     = abs(z);
                    a     = y;
                    a(y<yr)   = yr;
                    ya    = z./a;
                else %p = ya (phase map)

                    % IFFT
                    %             if isdataset(x)
                    %                 y   = real(ifft(x.data.*p));
                    %             else
                    y   = real(ifft(x.*p));   %assume x is class double
                    %            end
                    ya    = [];
                    p     = [];


                end
            end

            function [sp,imp]=pure(d,nr,f)
                % [sp,imp]=pure(d,nr,f)
                % sp purest row/column profiles
                % imp indexes of purest variables
                % d data matrix; nr (rank) number of pure components to search
                % if d(nspectra,nwave) imp gives purest nwave => sp are conc. profiles (nr,nspectra)
                % if d(nwave,nspectra) imp gives purest nspectra => sp are spectra profiles (nr,nwave)
                % f percent of noise allowed respect maximum of the average spectrum given in % (i.e. 1% or 0.1%))

                [nrow,ncol]=size(d);


                % calculation of the purity spectrum

                f=f/100;
                s=std(d);
                m=mean(d);
                ll=s.*s+m.*m;
                f=max(m)*f;
                p=s./(m+f);

                [mp,imp(1)]=max(p);

                % calculation of the correlation matrix
                % l=sqrt(m.*m+(s+f).*(s+f));

                l=sqrt((s.*s+(m+f).*(m+f)));
                % dl=d./(l'*ones(1,ncol));
                for j=1:ncol,
                    dl(:,j)=d(:,j)./l(j);
                end
                c=(dl'*dl)./nrow;

                % calculation of the weights
                % first weight

                w(1,:)=ll./(l.*l);
                p(1,:)=w(1,:).*p(1,:);
                s(1,:)=w(1,:).*s(1,:);
                % figure(1)
                % subplot(3,1,1),plot(m)
                % title('unweigthed mean, std and first pure spectrum')
                % subplot(3,1,2),plot(s)
                % subplot(3,1,3),plot(p)

                % pause

                % next weights


                for i=2:nr,

                    for j=1:ncol,
                        [dm]=wmat(c,imp,i,j);
                        w(i,j)=det(dm);
                        p(i,j)=p(1,j).*w(i,j);
                        s(i,j)=s(1,j).*w(i,j);
                    end

                    % next purest and standard deviation spectrum

                    % plot(p(i,:))
                    % plot(s(i,:))
                    % title('sd and purest spectrum')
                    % figure(i)
                    % subplot(2,1,1),plot(p(i,:))
                    % title('next pure spectrum and std dev. spectrum')
                    % subplot(2,1,2),plot(s(i,:))
                    % pause


                    [mp(i),imp(i)]=max(p(i,:));
                end


                for i=1:nr,
                    impi=imp(i);
                    sp(1:nrow,i)=d(1:nrow,impi);
                end

                % figure(nr+1)
                sp=normv2(sp');
                % plot(sp')

            end

            function [dm]=wmat(c,imp,irank,jvar)
                dm(1,1)=c(jvar,jvar);
                for k=2:irank,
                    kvar=imp(k-1);
                    dm(1,k)=c(jvar,kvar);
                    dm(k,1)=c(kvar,jvar);
                    for kk=2:irank,
                        kkvar=imp(kk-1);
                        dm(k,kk)=c(kvar,kkvar);
                    end
                end
            end

            function [sn]=normv2(s)
                % normalitzacio s=s/sqrt(sum(si)2))
                [m,n]=size(s);
                for i=1:m,
                    sr=sqrt(sum(s(i,:).*s(i,:)));
                    sn(i,:)=s(i,:)./sr;
                end
            end


            % ****************************** Subroutine****************************

            function [X] = makeXfromABC(A,B,C)

                na = size(A,1);
                nb = size(B,1);
                nc = size(C,1);
                ncomp = size(C,2);

                X = zeros(na,nb,nc);
                for i = 1:nc
                    Di = zeros(ncomp,ncomp);
                    for ii = 1:ncomp
                        Di(ii,ii) = C(i,ii);
                    end
                    X(:,:,i) = A*Di*B(:,:,i)';
                end

                X = permute(X,[3,2,1]);
            end

            function x=LSI(E,f,G,h);

                % x=LSI(E,f,G,h);
                % Solves problem LSI min(norm(Ex-f))) subject to Gx>h
                %
                % From Lawson & Hanson 74, p. 167
                %
                % Copyright 1998
                % Rasmus Bro
                % rasmus@optimax.dk/rb@kvl.dk

                [I,J]=size(E);
                r=rank(E);
                [u,s,v]=svd(E,0);

                K1=v(:,1:r);
                Q1=u(:,1:r);
                R=s(1:r,1:r);
                f1=Q1'*f;
                GKR=G*K1*pinv(R);
                z=ldp(GKR,h-GKR*f1);

                %z=Ry-f1;

                y=pinv(R)*(z+f1);
                x=K1*y;
                function b=ldp(G,h);

                    % find min||b|| subject to Gb=>h
                    %
                    % From Lawson & Hanson 74, p. 165
                    %
                    % Copyright 1998
                    % Rasmus Bro
                    % rasmus@optimax.dk/rb@kvl.dk


                    [I,J]=size(G);
                    E=[G';h'];
                    f=[zeros(J,1);1];
                    uhat=fastnnls(E'*E,E'*f);
                    r=E*uhat-f;

                    if norm(r)==0
                        disp(' No solution to LDP problem')
                    elseif r(J+1)~=0
                        b=-r(1:J)/r(J+1);
                    else
                        b=NaN;
                    end
                end
            end

            function [K, Pset] = fcnnls(C, A, CtC, CtA)
                % NNLS using normal equations and the fast combinatorial strategy
                %
                % I/O: [K, Pset] = fcnnls(C, A);
                % K = fcnnls(C, A);
                %
                % C is the nObs x lVar coefficient matrix
                % A is the nObs x pRHS matrix of observations
                % K is the lVar x pRHS solution matrix
                % Pset is the lVar x pRHS passive set logical array
                %
                % Pset: set of passive sets, one for each column
                % Fset: set of column indices for solutions that have not yet converged
                % Hset: set of column indices for currently infeasible solutions
                % Jset: working set of column indices for currently optimal solutions
                %
                % Implementation is based on [1] with bugfixes, direct passing of sufficient stats,
                % and preserving the active set over function calls.
                %
                % [1] Van Benthem, M. H., & Keenan, M. R. (2004). Fast algorithm for the
                %   solution of large‐scale non‐negativity‐constrained least squares problems.
                %   Journal of Chemometrics: A Journal of the Chemometrics Society, 18(10), 441-450.


                % Check the input arguments for consistency and initialize
                if nargin == 2
                    error(nargchk(2,2,nargin))
                    [nObs, lVar] = size(C);

                    if size(A,1)~= nObs, error('C and A have imcompatible sizes'), end
                    if size(C,1) == size(C,2)
                        %         warning('A square matrix "C" was input, ensure this is on purpose.')
                    end
                    pRHS = size(A,2);
                    % Precompute parts of pseudoinverse
                    CtC = C'*C; CtA = C'*A;
                else
                    [lVar,pRHS] = size(CtA);

                end

                if nargin == 2 && size(C,1) == size(C,2)
                    warning('fcnnls - The coefficient matrix (C) was square - is this true or are you passing C''C?')
                end

                W = zeros(lVar, pRHS);
                iter = 0;
                maxiter = 6*lVar;


                % Obtain the initial feasible solution and corresponding passive set
                K = cssls(CtC, CtA);
                Pset=K>0;
                K(~Pset) = 0;
                D=K;
                Fset = find(~all(Pset));

                % Active set algorithm for NNLS main loop
                iter_outer = 1;
                while ~isempty(Fset) && iter_outer < maxiter
                    iter_outer = iter_outer + 1;
                    % Solve for the passive variables (uses subroutine below)
                    K(:,Fset) = cssls(CtC, CtA(:,Fset), Pset(:,Fset));
                    % Find any infeasible solutions
                    %     Hset = Fset(find(any(K(:,Fset) < 0)));
                    Hset = Fset((any(K(:,Fset) < 0)));

                    % Make infeasible solutions feasible (standard NNLS inner loop)
                    if ~isempty(Hset)
                        nHset = length(Hset);
                        alpha = zeros(lVar, nHset);

                        while ~isempty(Hset) && (iter < maxiter)
                            iter = iter + 1;
                            alpha(:,1:nHset) = Inf;
                            % Find indices of negative variables in passive set
                            [i, j] = find(Pset(:,Hset) & (K(:,Hset) < 0));
                            hIdx = sub2ind([lVar nHset], i, j);
                            %             if length(i) ~= length(j)
                            %                 keyboard
                            %             end
                            %             negIdx = sub2ind(size(K), i, Hset(j)'); % org
                            negIdx = sub2ind(size(K), i, reshape(Hset(j),size(i)));  % jlh mod
                            alpha(hIdx) = D(negIdx)./(D(negIdx) - K(negIdx));
                            [alphaMin,minIdx] = min(alpha(:,1:nHset));
                            alpha(:,1:nHset) = repmat(alphaMin, lVar, 1);
                            D(:,Hset) = D(:,Hset)-alpha(:,1:nHset).*(D(:,Hset)-K(:,Hset));
                            idx2zero = sub2ind(size(D), minIdx, Hset);
                            D(idx2zero) = 0;
                            Pset(idx2zero) = 0;
                            K(:, Hset) = cssls(CtC, CtA(:,Hset), Pset(:,Hset));
                            Hset = find(any(K < 0));
                            nHset = length(Hset);
                        end
                    end
                    % Make sure the solution has converged
                    %if iter == maxiter, warning('Maximum number iterations exceeded'), end
                    % Check solutions for optimality
                    W(:,Fset) = CtA(:,Fset)-CtC*K(:,Fset);
                    Jset = find(all(~Pset(:,Fset).*W(:,Fset) <= 0));
                    Fset = setdiff(Fset, Fset(Jset));
                    % For non-optimal solutions, add the appropriate variable to Pset
                    if ~isempty(Fset)
                        [mx, mxidx] = max(~Pset(:,Fset).*W(:,Fset));
                        Pset(sub2ind([lVar pRHS], mxidx, Fset)) = 1;
                        D(:,Fset) = K(:,Fset);
                    end
                end

                % ****************************** Subroutine****************************
                function [K] = cssls(CtC, CtA, Pset)
                    % Solve the set of equations CtA = CtC*K for the variables in set Pset
                    % using the fast combinatorial approach
                    K = zeros(size(CtA));
                    if (nargin == 2) || isempty(Pset) || all(Pset(:))
                        K = CtC\CtA; % Not advisable if matrix is close to singular or badly scaled
                        %     K = pinv(CtC)*CtA;
                    else
                        [lVar, pRHS] = size(Pset);
                        codedPset = 2.^(lVar-1:-1:0)*Pset;
                        [sortedPset, sortedEset] = sort(codedPset);
                        breaks = diff(sortedPset);
                        breakIdx = [0 find(breaks) pRHS];
                        for k = 1:length(breakIdx)-1
                            cols2solve = sortedEset(breakIdx(k)+1:breakIdx(k+1));
                            vars = Pset(:,sortedEset(breakIdx(k)+1));
                            K(vars,cols2solve) = CtC(vars,vars)\CtA(vars,cols2solve);
                            %         K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
                        end
                    end
                end
            end
        end

        function [pks,locs_y,locs_x]=peaks2(data,varargin)
            % Find local peaks in 2D data.
            % Syntax chosen to be as close as possible to the original Matlab
            % 'findpeaks' function but not require any additional toolbox.
            %
            % SYNTAX:
            % pks=peaks(data) finds local peaks.
            %
            % [pks,locs_y,locs_x]=peaks(data) finds local peaks and their array
            % coordinates.
            %
            % [pks,locs_y,locs_x]=peaks(...,'MinPeakHeight',{scalar value}) only retains
            % those peaks which are equal to or greater than this absolute value.
            %
            % [pks,locs_y,locs_x]=peaks(...,'Threshold',{scalar value}) only retains
            % those peaks that are higher than their immediate surroundings by this value.
            %
            % [pks,locs_y,locs_x]=peaks(...,'MinPeakDistance',{scalar value}) finds peaks
            % separated by more than the specified minimum CARTESIAN peak distance (a
            % circle around the peak). It starts from the strongest peak and goes
            % iteratively lower. Any peak 'shadowed' in the vicinity of a stronger
            % peak is discarded.
            %
            % ALGORITHM:
            % A peak is considered to be a data point strictly greater than its
            % immediate neighbors. You can change this condition to 'greater or equal'
            % in the code, but be aware that in such case, it might create false
            % detections in flat areas, but these can be accounted for by introduction
            % of a small Threshold value.
            %
            % Even though this function is shared here free for use and any
            % modifications you might find useful, I would appreciate if you would
            % quote me in case you are going to use this function for any non-personal
            % tasks.
            % (C) Kristupas Tikuisis 2023.


            %% Initial data check
            % Let's simplify the function. Let it work on 1D or 2D data only, and check
            % if the input data satisfies this criteria:
            if ~ismatrix(data)
                error('Only 1D (vectors) or 2D (matrices) data accepted.');
            end


            %% Locate all peaks
            % A peak is a data point HIGHER than its immediate neighbors. There are 8
            % around eaxh point, and we will go through each. Oh yes, no escaping that.
            %
            % To introduce as little of intermediate variables and keep their memory
            % footprint as low as possible, I will introduce 2 logical arrays: one to
            % mark the peaks and be iteratively updated until we check all its
            % neighbors; and another temporal variable just to prepare data for
            % comparison (mind about the edge points which do not have any neighbors!):
            ispeak=false(size(data)); % to store peak flags.
            isgreater=true(size(data)); % to store comparison result for one particular neighbor.

            % Now start analyzing every data point.
            %
            % 1st neighbor immediatelly to the left:
            ispeak=([true(size(data,1),1) [data(:,2:end)>data(:,1:end-1)]]); % for this case, we can update the peak array directly.
            %
            % 2nd neighbor at the top-left:
            isgreater(2:end,2:end)=(data(2:end,2:end)>data(1:end-1,1:end-1)); % this time, due to points on the diagonal, a temporary array will have to be involved.
            ispeak=ispeak&isgreater; % time to update the peak array.
            %
            % 3rd neighbor immediatelly at the top:
            ispeak=ispeak&([true(1,size(data,2)); (data(2:end,:)>data(1:end-1,:))]); % once again, for this case, we can update the peak array directly.
            %
            % 4th neighbor at the top-right:
            isgreater=true(size(data)); % rebuild a fresh array for comparison...
            isgreater(2:end,1:end-1)=(data(2:end,1:end-1)>data(1:end-1,2:end));
            ispeak=ispeak&isgreater;
            %
            % 5th neighbor immediatelly to the right:
            ispeak=ispeak&([(data(:,1:end-1)>data(:,2:end)) true(size(data,1),1)]); % once again, for this case, we can update the peak array directly.
            %
            % 6th neighbor to the bottom-right:
            isgreater=true(size(data));
            isgreater(1:end-1,1:end-1)=(data(1:end-1,1:end-1)>data(2:end,2:end));
            ispeak=ispeak&isgreater;
            %
            % 7th neighbor immediatelly at the bottom:
            ispeak=ispeak&([(data(1:end-1,:)>data(2:end,:)); true(1,size(data,2))]); % once again, for this case, we can update the peak array directly.
            %
            % 8th neighbor at the bottom-left:
            isgreater=true(size(data));
            isgreater(1:end-1,2:end)=(data(1:end-1,2:end)>data(2:end,1:end-1));
            ispeak=ispeak&isgreater;
            %
            % Discard the temporary variable:
            clear isgreater
            % By now, raw peak indentification is completed.


            %% Return final results
            % First, essential raw peak location steps. Peak values and locations:
            locs=find(ispeak); %clear ispeak
            pks=data(locs);
            % Note that LINEAR indices have been returned. Let's turn them to array
            % indices for final output:
            [locs_y,locs_x]=ind2sub(size(data),locs);


            %% Perform customary post-processing determined by optional function parameters.
            % First, check if any parameters were supplied:
            if isempty(varargin)
                return
            end
            % ...then a quality check - the parameters should come in pairs, therefore
            % the length should be even:
            if mod(length(varargin),2)~=0
                warning('Optional name-value parameters should come in pairs. Something is missing. Peak search will proceed with default values.');
                return
            end
            % ...after this step, we can sort the input into names (parameters) and
            % values:
            params=varargin(1:2:end);
            vals=varargin(2:2:end);
            % ...last quality check - all even values should be CHAR entries specifying
            % a parameter to be adjusted:
            ischarparam=cellfun(@ischar,params);
            if any(~ischarparam)
                warning('Some parameters are not written as characters and not recognisable. They should always come is name(char)-value pairs. Peak search will proceed with default values.');
                return
            end; clear ischarparam varargin

            %--------------------------------------------------------------------------
            % No go over all supplied parameters and check the found peaks accordingly.

            % 1. MinPeakHeight - absolute minimum value for a peak:
            isparm=find(cellfun(@(x)isequal(x,'MinPeakHeight'),params),1);
            if ~isempty(isparm)
                % Locate which peaks satisfy this condition:
                suitable=(pks>=vals{isparm});
                % ...and only keep those:
                pks=pks(suitable);
                locs=locs(suitable);
                clear suitable
            end

            % 2. Threshold - peak must be greater than its neighbours by this value.
            isparm=find(cellfun(@(x)isequal(x,'Threshold'),params),1);
            if ~isempty(isparm)
                % For this, we will need to convert linear indices to array indices:
                [row_y,col_x]=ind2sub(size(data),locs);
                % These will be the original indices.

                % This is how array coordinates would change relatively around each
                % peak (y,x):
                % (-1,-1)  (-1,0)  (-1,+1)
                % ( 0,-1)  ( 0,0)  ( 0,+1)
                % (+1,-1)  (+1,0)  (+1,+1)
                % ...turned into vectors disregarding the (0,0), the center data point:
                delta_y=[-1 -1 -1 0 0 +1 +1 +1];
                delta_x=[-1 0 +1 -1 +1 -1 0 +1];
                % Let's add these deltas to the detected peak positions to get the
                % coordinates of their immediate neighbors:
                neighbor_locs_y=row_y+delta_y;
                neighbor_locs_x=col_x+delta_x; clear row_y col_x delta_x delta_y
                % ...don't forget to check for unrealistic indices beyond array
                % borders:
                neighbor_locs_y(neighbor_locs_y<1)=1;
                neighbor_locs_y(neighbor_locs_y>size(data,1))=size(data,1);
                neighbor_locs_x(neighbor_locs_x<1)=1;
                neighbor_locs_x(neighbor_locs_x>size(data,2))=size(data,2);
                % ...convert to linear indices:
                neighbor_locs=sub2ind(size(data),neighbor_locs_y,neighbor_locs_x);
                clear neighbor_locs_y neighbor_locs_x

                % So we have neighbor values by now. Are our peaks higher than those by
                % the set Threshold value?
                suitable=(data(locs)-vals{isparm}>=data(neighbor_locs));
                % Now check for those cases when by mistake (earlier step for checking
                % for indices beyond array boundaries) a peak itself is taken as a
                % neighbor as well:
                suitable(data(neighbor_locs)==data(locs))=true;

                % Only keep those is they are greater than ALL neighbors (in other
                % words, those elements where NONE are lesser):
                suitable=~any(~suitable,2); clear neighbor_locs

                % Final step - locate suitable element indices:
                suitable=find(suitable);

                % That's it, modify the output array:
                pks=pks(suitable);
                locs=locs(suitable); clear suitable
            end

            % 3. 'MinPeakDistance'
            isparm=find(cellfun(@(x)isequal(x,'MinPeakDistance'),params),1);
            if ~isempty(isparm)

                % First, sort the peaks in order of amplitude:
                [pks_sorted,idx]=sort(pks,'descend');
                locs_sorted=locs(idx); clear idx

                % The flow is as follows: start from the highest peak and discard any
                % other peaks closer than the CARTESIAN MinPeakDistance (that is, the
                % CIRCLE around the peak is going to be checked); then continue until
                % the whole list (updated iteratively as items might get removed) has
                % been checked.

                % Convert locations to array indices:
                [row_y,col_x]=ind2sub(size(data),locs_sorted);

                % Start from the highest peak:
                this_peak=1;
                while this_peak<(length(pks_sorted)+1)

                    % Cartesian distances to ALL its remaining & yet unchecked neighbors
                    % (including itself):
                    dist=sqrt((row_y-row_y(this_peak)).^2+(col_x-col_x(this_peak)).^2);
                    % Now simply check which neighbors are WITHIN the MinPeakDistance
                    % BUT also nonzero (the peak should not be compared to itself):
                    within=( (dist<=vals{isparm}) & (dist~=0) );
                    % ...and delete those entries satisfying the condition:
                    pks_sorted(within)=[];
                    locs_sorted(within)=[];
                    row_y(within)=[]; col_x(within)=[];

                    % Update the peak counter:
                    this_peak=this_peak+1;

                end; clear this_peak within dist row_y col_x

                % Update the peak location and value list:
                pks=pks_sorted;
                locs=locs_sorted;
            end

            % I haven't figured out how to do it more elegantly, but turn the default
            % linear indices to array indices:
            [locs_y,locs_x]=ind2sub(size(data),locs);

            %==========================================================================
            % End of the entire function
        end

        function [Fac11,Fac12,Fac13,Fac21,Fac22,Fac23] = get_rank_of_modes(outputstruc,filtered)
            if filtered == 0
                % do for mz1
                data = double(outputstruc.mz1_temp);

                mode1 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode2 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode3 = reshape(data,size(data,1),size(data,2)*size(data,3));

                s1 = svd(mode1,'econ');
                s2 = svd(mode2,'econ');
                s3 = svd(mode3,'econ');


                exVar1 = s1.^2./sum(s1.^2);
                exVar2 = s2.^2./sum(s2.^2);
                exVar3 = s3.^2./sum(s3.^2);

                i = 1;
                while sum(exVar1(1:i))/sum(exVar1) < 0.999
                    i = i+1;
                end
                Fac11 = i;

                i = 1;
                while sum(exVar2(1:i))/sum(exVar2) < 0.999
                    i = i+1;
                end
                Fac12 = i;

                i = 1;
                while sum(exVar3(1:i))/sum(exVar3) < 0.999
                    i = i+1;
                end
                Fac13 = i;

                % do for mz2
                data = double(outputstruc.mz2_temp);

                mode1 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode2 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode3 = reshape(data,size(data,1),size(data,2)*size(data,3));

                s1 = svd(mode1,'econ');
                s2 = svd(mode2,'econ');
                s3 = svd(mode3,'econ');


                exVar1 = s1.^2./sum(s1.^2);
                exVar2 = s2.^2./sum(s2.^2);
                exVar3 = s3.^2./sum(s3.^2);

                i = 1;
                while sum(exVar1(1:i))/sum(exVar1) < 0.999
                    i = i+1;
                end
                Fac21 = i;

                i = 1;
                while sum(exVar2(1:i))/sum(exVar2) < 0.999
                    i = i+1;
                end
                Fac22 = i;

                i = 1;
                while sum(exVar3(1:i))/sum(exVar3) < 0.999
                    i = i+1;
                end
                Fac23 = i;

                %%
            else
 
                % do for mz1
                data = double(outputstruc.mz1_filtered);

                mode1 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode2 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode3 = reshape(data,size(data,1),size(data,2)*size(data,3));

                s1 = svd(mode1,'econ');
                s2 = svd(mode2,'econ');
                s3 = svd(mode3,'econ');


                exVar1 = s1.^2./sum(s1.^2);
                exVar2 = s2.^2./sum(s2.^2);
                exVar3 = s3.^2./sum(s3.^2);

                i = 1;
                while sum(exVar1(1:i))/sum(exVar1) < 0.99
                    i = i+1;
                end
                Fac11 = i;

                i = 1;
                while sum(exVar2(1:i))/sum(exVar2) < 0.99
                    i = i+1;
                end
                Fac12 = i;

                i = 1;
                while sum(exVar3(1:i))/sum(exVar3) < 0.99
                    i = i+1;
                end
                Fac13 = i;

                % do for mz2
                data = double(outputstruc.mz2_filtered);

                mode1 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode2 = reshape(data,size(data,1),size(data,2)*size(data,3));
                data = permute(data,[2,3,1]);
                mode3 = reshape(data,size(data,1),size(data,2)*size(data,3));

                s1 = svd(mode1,'econ');
                s2 = svd(mode2,'econ');
                s3 = svd(mode3,'econ');


                exVar1 = s1.^2./sum(s1.^2);
                exVar2 = s2.^2./sum(s2.^2);
                exVar3 = s3.^2./sum(s3.^2);

                i = 1;
                while sum(exVar1(1:i))/sum(exVar1) < 0.99
                    i = i+1;
                end
                Fac21 = i;

                i = 1;
                while sum(exVar2(1:i))/sum(exVar2) < 0.99
                    i = i+1;
                end
                Fac22 = i;

                i = 1;
                while sum(exVar3(1:i))/sum(exVar3) < 0.99
                    i = i+1;
                end
                Fac23 = i;
               
            end

        end

        function [sp,Xhat] = createOutputFromModel(model,ind_target_mz)

            if iscell(model)
                dim = max(size(model));
                for i = 1:dim
                    fits(i) = model{i}.detail.fit.PercVar;
                end
                [val,ind1] = sort(fits,'descend');
                model_selected = model{ind1(1)};
                [val,ind2] = max(model_selected.spectra(ind_target_mz,:));
                sp = model_selected.spectra(:,ind2);
                Xhat = makeXfromABC(model_selected.spectra(:,ind2),model_selected.elutionprofiles(:,ind2,:),model_selected.scores(:,ind2));
            end
            function [X] = makeXfromABC(A,B,C)

                na = size(A,1);
                nb = size(B,1);
                nc = size(C,1);
                ncomp = size(C,2);

                X = zeros(na,nb,nc);
                for i = 1:nc
                    Di = zeros(ncomp,ncomp);
                    for ii = 1:ncomp
                        Di(ii,ii) = C(i,ii);
                    end
                    X(:,:,i) = A*Di*B(:,:,i)';
                end

                X = permute(X,[3,2,1]);
            end
        end

        function make_result_tables(outputstruc,outputFile,mtype)
            mz_target = outputstruc.mz_target;
            if mtype == 1
                mz_axis1 = outputstruc.axis_mz1_filtered;
                mz_axis2 = outputstruc.axis_mz2_filtered;
            elseif mtype == 2
                mz_axis1 = outputstruc.mz1_axis_temp;
                mz_axis2 = outputstruc.mz2_axis_temp;
            end
            % Assume that you already have a vector 'mz_axis'
            if length(mz_axis1) > length(mz_axis2)
                mz1 = zeros(length(mz_axis1), 1); % Initialize column 'mz1' with zeros
                mz2 = zeros(length(mz_axis1), 1); % Initialize column 'mz2' with zeros
                int1 = zeros(length(mz_axis1), 1); % Initialize column 'mz1' with zeros
                int2 = zeros(length(mz_axis1), 1); % Initialize column 'mz1' with zeros
            else
                mz1 = zeros(length(mz_axis2), 1); % Initialize column 'mz1' with zeros
                mz2 = zeros(length(mz_axis2), 1); % Initialize column 'mz2' with zeros
                int1 = zeros(length(mz_axis2), 1); % Initialize column 'mz1' with zeros
                int2 = zeros(length(mz_axis2), 1); % Initialize column 'mz1' with zeros
            end

            keep_mz1 = mz_axis1 < mz_target+cutoff;
            keep_mz2 = mz_axis2 < mz_target+cutoff;

            int1(1:length(outputstruc.mz1_mass_spectrum(keep_mz1))) = outputstruc.mz1_mass_spectrum(keep_mz1);
            int1 = int1./(max(int1))*1000;

            int2(1:length(outputstruc.mz2_mass_spectrum(keep_mz2))) = outputstruc.mz2_mass_spectrum(keep_mz2);
            int2 = int2./(max(int2))*1000;

            mz1(1:length(mz_axis1(keep_mz1))) = mz_axis1(keep_mz1);
            mz2(1:length(mz_axis2(keep_mz2))) = mz_axis2(keep_mz2);


            [val,ind] = sort(int1,'descend');
            mz1 = mz1(ind);
            int1 = val;


            [val,ind] = sort(int2,'descend');
            mz2 = mz2(ind);
            int2 = val;

            % Create a table with columns 'mz1' and 'mz2
            result_table = table(mz1,int1,mz2,int2);

            % Export the table to an Excel file
            writetable(result_table, outputFile);
        end

        function PeakApexMethod(inputStruc,outputFile)

            %             [~,xmz1] = max(inputStruc.mz1_el_profil_lc2);
            %             [~,ymz1] = max(inputStruc.mz1_el_profil_lc1);
            %
            %             [~,xmz2] = max(inputStruc.mz2_el_profil_lc2);
            %             [~,ymz2] = max(inputStruc.mz2_el_profil_lc1);

            [~,ind] = min(abs(inputStruc.mz_target-inputStruc.axis_mz1_temp));
            [ints,y1,x1] = peaks2(squeeze(double(inputStruc.mz1_temp(:,:,ind))));
            [~,ind] = sort(ints,'descend');
            y1 = y1(ind(1));
            x1 = x1(ind(1));

            [~,ind] = min(abs(inputStruc.mz_target-inputStruc.axis_mz2_temp));
            [ints,y2,x2] = peaks2(squeeze(double(inputStruc.mz2_temp(:,:,ind))));
            [~,ind] = sort(ints,'descend');
            y2 = y2(ind(1));
            x2 = x2(ind(1));


            int1 = squeeze(double(inputStruc.mz1_temp(y1,x1,:)));
            int2 = squeeze(double(inputStruc.mz2_temp(y2,x2,:)));

            mz1  = double(inputStruc.axis_mz1_temp(int1>max(int1)*0.001));
            mz2  = double(inputStruc.axis_mz2_temp(int2>max(int2)*0.001));

            int1 = int1(int1>max(int1)*0.001);
            int2 = int2(int2>max(int2)*0.001);

            int1 = int1(mz1 < (inputStruc.mz_target+cutoff));
            int2 = int2(mz2 < (inputStruc.mz_target+cutoff));

            int1 = int1./max(int1)*1000;
            int2 = int2./max(int2)*1000;

            mz1  = mz1(mz1 < (inputStruc.mz_target+cutoff));
            mz2  = mz2(mz2 < (inputStruc.mz_target+cutoff));


            [~,ind] = sort(int1,'descend');
            int1 = int1(ind);
            mz1 = mz1(ind);

            [~,ind] = sort(int2,'descend');
            int2 = int2(ind);
            mz2 = mz2(ind);


            if length(int1) > length(int2)
                deltaint = length(int1)-length(int2);
                int2 = [int2; zeros(deltaint,1)];
                mz2  = [mz2; zeros(deltaint,1)];
            else
                deltaint = length(int2)-length(int1);
                int1 = [int1; zeros(deltaint,1)];
                mz1  = [mz1; zeros(deltaint,1)];
            end
            result_table = table(mz1,int1,mz2,int2);

            % Export the table to an Excel file
            writetable(result_table, outputFile)
        end

        function [cosineSimilarity,shared_fragments,all_shared_fragments] = refspectra_analysis(mz_threshold)
            if isfile("ref_spectrum.xlsx")
                ref_spectrum = readtable("ref_spectrum.xlsx",Sheet="RefSpec");
                mz_dif = mz_threshold;
                if ismember('SpecIndex', ref_spectrum.Properties.VariableNames)
                    if ~isnan(ref_spectrum.SpecIndex)
                        num_spectra =  unique(ref_spectrum.SpecIndex);
                    else

                        ref_spectrum.SpecIndex = ones(size(ref_spectrum,1),1);
                        num_spectra = 1;
                    end
                else
                    ref_spectrum.SpecIndex = ones(size(ref_spectrum,1),1);
                    num_spectra = 1;

                end

                for num = num_spectra'
                    clear ind Ref_spectrum_neu
                    all_i = [];
                    ind_logical = ref_spectrum.SpecIndex == num;
                    ref_spectrum_temp.m_z = ref_spectrum.m_z(ind_logical);
                    ref_spectrum_temp.Rel_Int = ref_spectrum.Rel_Int(ind_logical);


                    MF_spectrum = readtable("spectra_MF.xlsx",Sheet="Sheet1");
                    MF_spectrum.mz2(MF_spectrum.int2==0) = 0;
                    MF_SIT_spectrum = readtable("spectra_MF_SIT.xlsx",Sheet="Sheet1");
                    MF_SIT_spectrum.mz2(MF_SIT_spectrum.int2==0) = 0;
                    PAM_spectrum = readtable("spectra_PAM.xlsx",Sheet="Sheet1");
                    SIT_spectrum = readtable("spectra_SIT.xlsx",Sheet="Sheet1");
                    SIT_spectrum.mz2(SIT_spectrum.int2==0) = 0;
                    % PAM is the longest vector - set as standard
                    max_size = size(PAM_spectrum,1);

                    %adjust length of all other vectors to the same length and sort the
                    %elements to fit with axis position of PAM
                    MF_SIT_spectrum_neu = zeros(max_size,2);
                    [val,ind1] = intersect(round(PAM_spectrum.mz2,4),round(MF_SIT_spectrum.mz2,4));
                    [val,ind2] = intersect(round(MF_SIT_spectrum.mz2,4),val);

                    MF_SIT_spectrum_neu(ind1,1) = MF_SIT_spectrum.mz2(ind2);
                    MF_SIT_spectrum_neu(ind1,2) = MF_SIT_spectrum.int2(ind2);
                    clear val ind1 ind2

                    MF_spectrum_neu = zeros(max_size,2);
                    [val,ind1] = intersect(PAM_spectrum.mz2,MF_spectrum.mz2);
                    [val,ind2] = intersect(MF_spectrum.mz2,val);

                    MF_spectrum_neu(ind1,1) = MF_spectrum.mz2(ind2);
                    MF_spectrum_neu(ind1,2) = MF_spectrum.int2(ind2);
                    clear val ind1 ind2

                    PAM_spectrum_neu(:,1) = PAM_spectrum.mz2;
                    PAM_spectrum_neu(:,2) = PAM_spectrum.int2;


                    SIT_spectrum_neu = zeros(max_size,2);
                    [val,ind1] = intersect(PAM_spectrum.mz2,SIT_spectrum.mz2);
                    [val,ind2] = intersect(SIT_spectrum.mz2,val);

                    SIT_spectrum_neu(ind1,1) = SIT_spectrum.mz2(ind2);
                    SIT_spectrum_neu(ind1,2) = SIT_spectrum.int2(ind2);
                    clear val ind1 ind2

                    % intialize lenght of Ref_spectrum_neu
                    Ref_spectrum_neu = zeros(max_size,2);



                    % PAM spectrum and preparation of reference spectrum
                    iii = 0;
                    for ii = 1:size(ref_spectrum_temp.m_z,1)

                        [val,ind] = min(abs(ref_spectrum_temp.m_z(ii)-PAM_spectrum_neu(:,1)));
                        if val <= mz_dif
                            iii = iii+1;
                            ref_spectrum_temp.m_z(iii) = PAM_spectrum_neu(ind,1);
                            ref_spectrum_temp.Rel_Int(iii)  = ref_spectrum_temp.Rel_Int(ii);
                            all_i(iii) = ind;
                        end
                    end
                    clear val ind
                    [val,ind] = intersect(ref_spectrum_temp.m_z,PAM_spectrum_neu(unique(all_i),1));
                    for iiii = 1:length(ind)
                        ind_temp = [];
                        ind_temp = find(ref_spectrum_temp.m_z(ind(iiii)) == PAM_spectrum_neu(:,1));
                        Ref_spectrum_neu(ind_temp,1) = ref_spectrum_temp.m_z(ind(iiii));
                        Ref_spectrum_neu(ind_temp,2) = ref_spectrum_temp.Rel_Int(ind(iiii));
                    end
                    Ref_spectrum_neu(:,2) = Ref_spectrum_neu(:,2)./norm(Ref_spectrum_neu(:,2),'fro');
                    PAM_spectrum_neu(:,2) = PAM_spectrum_neu(:,2)./norm(PAM_spectrum_neu(:,2),'fro');

                    %                     cosineSimilarity(num,1) = Ref_spectrum_neu(:,2)'*PAM_spectrum_neu(:,2);
                    if exist("all_i")
                        cosineSimilarity(num,1) = Ref_spectrum_neu(:,2)'*PAM_spectrum_neu(:,2);
                        shared_fragments(num,1) = length(unique(all_i));
                        all_shared_fragments{num,1} = all_i;
                    else
                        cosineSimilarity(num,1) = 0;
                        shared_fragments(num,1) = 0;
                        all_shared_fragments{num,1} = 0;
                    end
                    %                     f = figure;
                    %                     subplot(2,2,1)
                    %                     stem(PAM_spectrum_neu(:,1),PAM_spectrum_neu(:,2),'Marker','none','Color',[0.6,0.6,0.6],'LineWidth',2)
                    %                     hold on;
                    %                     if exist("all_i")
                    %                         stem(PAM_spectrum_neu(all_i,1),PAM_spectrum_neu(all_i,2),'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2)
                    %                         hold on;
                    %                     end
                    %                     stem(Ref_spectrum_neu(:,1),Ref_spectrum_neu(:,2)*-1,'Marker','none','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
                    %                     ylim([-1,1])
                    clear ind val all_i ind_temp
                    % SIT spectrum
                    iii = 0;
                    for ii = 1:size(ref_spectrum_temp.m_z,1)

                        [val,ind] = min(abs(ref_spectrum_temp.m_z(ii)-SIT_spectrum_neu(:,1)));
                        if val <= mz_dif
                            iii = iii+1;
                            all_i(iii) = ind;
                        end
                    end


                    SIT_spectrum_neu(:,2) = SIT_spectrum_neu(:,2)./norm(SIT_spectrum_neu(:,2),'fro');

                    %                     cosineSimilarity(num,2) = Ref_spectrum_neu(:,2)'*SIT_spectrum_neu(:,2);
                    if exist("all_i")
                        cosineSimilarity(num,2) = Ref_spectrum_neu(:,2)'*SIT_spectrum_neu(:,2);
                        shared_fragments(num,2) = length(unique(all_i));
                        all_shared_fragments{num,2} = all_i;
                    else
                        all_shared_fragments{num,2} = 0;
                        cosineSimilarity(num,2) = 0;
                        shared_fragments(num,2) = 0;
                    end
                    %                     subplot(2,2,2)
                    %                     stem(SIT_spectrum_neu(:,1),SIT_spectrum_neu(:,2),'Marker','none','Color',[0.6,0.6,0.6],'LineWidth',2)
                    %                     hold on;
                    %                     if exist("all_i")
                    %                         stem(SIT_spectrum_neu(all_i,1),SIT_spectrum_neu(all_i,2),'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2)
                    %                     end
                    %                     hold on;
                    %                     stem(Ref_spectrum_neu(:,1),Ref_spectrum_neu(:,2)*-1,'Marker','none','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
                    %                     ylim([-1,1])
                    clear val ind all_i

                    % MF spectrum
                    iii = 0;
                    for ii = 1:size(ref_spectrum_temp.m_z,1)
                        [val,ind] = min(abs(ref_spectrum_temp.m_z(ii)-MF_spectrum_neu(:,1)));

                        if val <= mz_dif
                            iii = iii+1;
                            all_i(iii) = ind;
                        end
                    end

                    MF_spectrum_neu(:,2) = MF_spectrum_neu(:,2)./norm(MF_spectrum_neu(:,2),'fro');

                    %                     cosineSimilarity(num,3) = Ref_spectrum_neu(:,2)'*MF_spectrum_neu(:,2);
                    if exist("all_i")
                        cosineSimilarity(num,3) = Ref_spectrum_neu(:,2)'*MF_spectrum_neu(:,2);
                        shared_fragments(num,3) = length(unique(all_i));
                        all_shared_fragments{num,3} = all_i;
                    else
                        all_shared_fragments{num,3} = 0;
                        cosineSimilarity(num,3) = 0;
                        shared_fragments(num,3) = 0;
                    end
                    %                     subplot(2,2,3)
                    %                     stem(MF_spectrum_neu(:,1),MF_spectrum_neu(:,2),'Marker','none','Color',[0.6,0.6,0.6],'LineWidth',2)
                    %                     hold on;
                    %                     if exist("all_i")
                    %                         stem(MF_spectrum_neu(all_i,1),MF_spectrum_neu(all_i,2),'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2)
                    %                     end
                    %                     hold on;
                    %                     stem(Ref_spectrum_neu(:,1),Ref_spectrum_neu(:,2)*-1,'Marker','none','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
                    %                     ylim([-1,1])
                    clear val ind all_i

                    % MF and SIT together
                    iii = 0;
                    for ii = 1:size(ref_spectrum_temp.m_z,1)
                        [val,ind] = min(abs(ref_spectrum_temp.m_z(ii)-MF_SIT_spectrum_neu(:,1)));
                        if val <= mz_dif
                            iii = iii+1;
                            all_i(iii) = ind;
                        end
                    end


                    MF_SIT_spectrum_neu(:,2) = MF_SIT_spectrum_neu(:,2)./norm(MF_SIT_spectrum_neu(:,2),'fro');

                    %                     cosineSimilarity(num,4) = Ref_spectrum_neu(:,2)'*MF_SIT_spectrum_neu(:,2);
                    if exist("all_i")
                        cosineSimilarity(num,4) = Ref_spectrum_neu(:,2)'*MF_SIT_spectrum_neu(:,2);
                        shared_fragments(num,4) = length(unique(all_i));
                        all_shared_fragments{num,4} = all_i;
                    else
                        cosineSimilarity(num,4) = 0;
                        shared_fragments(num,4) = 0;
                        all_shared_fragments{num,4} = 0;
                    end
                    %                     subplot(2,2,4)
                    %                     stem(MF_SIT_spectrum_neu(:,1),MF_SIT_spectrum_neu(:,2),'Marker','none','Color',[0.6,0.6,0.6],'LineWidth',2)
                    %                     hold on;
                    %                     if exist("all_i")
                    %                         stem(MF_SIT_spectrum_neu(all_i,1),MF_SIT_spectrum_neu(all_i,2),'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2)
                    %                     end
                    %                     hold on;
                    %                     stem(Ref_spectrum_neu(:,1),Ref_spectrum_neu(:,2)*-1,'Marker','none','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
                    %                     ylim([-1,1])
                    %                     saveas(f, ['compareWithref' num2str(num) '.svg'], 'svg');
                    %
                    %                     close all
                    clear val ind all_i


                end
            end


        end

    end

end

