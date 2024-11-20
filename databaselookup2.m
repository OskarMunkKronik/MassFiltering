
sus = readtable('Data\SuspectList.xlsx');

for ii = 1:size(sus,1)

database_name = "Data\MSMS-Public_experimentspectra-pos-VS19.msp";
excel_file = strcat("ref_spectra_",num2str(ii),"_.XLSX");
compound_name = sus.Name{ii};
try
extract_spectra_from_msp(database_name,compound_name,excel_file);
end

end

function extract_spectra_from_msp(msp_file, compound_name, output_excel_file)
    % Open the MSP file for reading
    fid = fopen(msp_file, 'r');
    if fid == -1
        error('Could not open the file: %s', msp_file);
    end
    
    % Initialize variables
    isCompoundFound = false; % Flag to track when the compound name is found
    spectra = []; % To store spectra data
    all_spectra = []; % To store all spectra for output
    spec_index = 0; % Index to track spectrum number
    
    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid); % Get the current line
        
        % Check if the line starts with 'NAME:' and matches the compound name
        if startsWith(line, 'NAME:')
            % Extract the compound name from the line
            current_name = strtrim(strrep(line, 'NAME:', ''));
            
            if strcmpi(current_name, compound_name)
                isCompoundFound = true;
                spectra = []; % Reset spectra for this entry
                spec_index = spec_index + 1; % Increment spectrum index
            else
                isCompoundFound = false; % Reset if the name does not match
            end
        end
        
        % Check for 'Num Peaks' line
        if isCompoundFound && startsWith(line, 'Num Peaks:')
            % Get the number of peaks (not strictly needed here)
            num_peaks = str2double(strrep(line, 'Num Peaks:', ''));
        end
        
        % If peaks data is found (pairs of m/z and intensity), parse it
        if isCompoundFound && ~isempty(line) && ~startsWith(line, 'NAME:') && ~startsWith(line, 'Num Peaks:')
            % Try to parse the m/z and intensity values
            data = str2num(line); %#ok<ST2NM>
            if ~isempty(data) && numel(data) == 2
                % Append the m/z, intensity, and spectrum index to the output
                all_spectra = [all_spectra; data(1), data(2), spec_index]; %#ok<AGROW>
            end
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Check if any spectra were found
    if isempty(all_spectra)
        error('No spectra found for the compound: %s', compound_name);
    end
    
    % Define the headers for the Excel output
    headers = {'m/z', 'Rel_Int', 'SpecIndex'};
    T = array2table(all_spectra,'VariableNames',headers);
    % Write the headers and the data to the output Excel file
    writetable(T,output_excel_file,"FileType","spreadsheet","Sheet","RefSpec")
    
    fprintf('Mass spectra extracted and saved to %s\n', output_excel_file);
end
