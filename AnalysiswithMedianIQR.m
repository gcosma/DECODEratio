clear
% Specify the folder paths
baseFolder = 'C:\Users\cogc\OneDrive - Loughborough University\MATLAB\MATLAB\ALLDataDatesHere\PERFECTCODE\SAIL\';

%change 3 things - 1
analysisFolderPath = fullfile(baseFolder, 'males');%remember to chage the condition number
%analysisFolderPath = fullfile(baseFolder, 'females');
dataFolder = baseFolder;

% Create the analysis folder if it doesn't exist
if ~exist(analysisFolderPath, 'dir')
    mkdir(analysisFolderPath);
    disp(['Created analysis folder: ', analysisFolderPath]);
else
    disp(['Analysis folder already exists: ', analysisFolderPath]);
end

%2
% Load the data
%load(fullfile(dataFolder, 'femaledatebirth.mat'), 'datamerge');
load(fullfile(dataFolder, 'maledatebirth.mat'), 'datamerge');

%3
% Get condition names (columns 1 to 40)
%conditions = datamerge.Properties.VariableNames(1:40);%females
conditions = datamerge.Properties.VariableNames(1:38); %males

% Initialize cell arrays for storing results
diagnosesBelow45 = cell(height(datamerge), length(conditions));
diagnoses45to64 = cell(height(datamerge), length(conditions));
diagnoses65AndAbove = cell(height(datamerge), length(conditions));

% Process each patient
for i = 1:height(datamerge)
    yearOfBirth = datamerge{i, end};
    
    for j = 1:length(conditions)
        diagnosisDate = datamerge{i, j}{1};  % Access the string inside the cell
        if ~isempty(diagnosisDate) && ~isempty(strtrim(diagnosisDate))
            diagDate = datetime(diagnosisDate, 'InputFormat', 'dd/MM/yyyy');
            ageAtDiagnosis = year(diagDate) - yearOfBirth;
            if ageAtDiagnosis < 45
                diagnosesBelow45{i, j} = diagnosisDate;
                diagnoses45to64{i, j} = '';
                diagnoses65AndAbove{i, j} = '';
            elseif ageAtDiagnosis >= 45 && ageAtDiagnosis < 65
                diagnosesBelow45{i, j} = '';
                diagnoses45to64{i, j} = diagnosisDate;
                diagnoses65AndAbove{i, j} = '';
            else
                diagnosesBelow45{i, j} = '';
                diagnoses45to64{i, j} = '';
                diagnoses65AndAbove{i, j} = diagnosisDate;
            end
        else
            diagnosesBelow45{i, j} = '';
            diagnoses45to64{i, j} = '';
            diagnoses65AndAbove{i, j} = '';
        end
    end
end

% Convert cell arrays to tables
diagnosesBelow45 = cell2table(diagnosesBelow45, 'VariableNames', conditions);
diagnoses45to64 = cell2table(diagnoses45to64, 'VariableNames', conditions);
diagnoses65AndAbove = cell2table(diagnoses65AndAbove, 'VariableNames', conditions);

% Add a column for Year of Birth
diagnosesBelow45.YearOfBirth = datamerge{:, end};
diagnoses45to64.YearOfBirth = datamerge{:, end};
diagnoses65AndAbove.YearOfBirth = datamerge{:, end};

% Save the results
save(fullfile(analysisFolderPath, 'diagnoses_below_45.mat'), 'diagnosesBelow45');
save(fullfile(analysisFolderPath, 'diagnoses_45_to_64.mat'), 'diagnoses45to64');
save(fullfile(analysisFolderPath, 'diagnoses_65_and_above.mat'), 'diagnoses65AndAbove');
writetable(diagnosesBelow45, fullfile(analysisFolderPath, 'diagnoses_below_45.csv'));
writetable(diagnoses45to64, fullfile(analysisFolderPath, 'diagnoses_45_to_64.csv'));
writetable(diagnoses65AndAbove, fullfile(analysisFolderPath, 'diagnoses_65_and_above.csv'));

% Display summary
totalPatients = height(diagnosesBelow45);
patientsBelow45 = sum(any(~cellfun(@isempty, table2cell(diagnosesBelow45(:, 1:end-1))), 2));
patients45to64 = sum(any(~cellfun(@isempty, table2cell(diagnoses45to64(:, 1:end-1))), 2));
patients65AndAbove = sum(any(~cellfun(@isempty, table2cell(diagnoses65AndAbove(:, 1:end-1))), 2));
disp(['Total number of patients: ' num2str(totalPatients)]);
disp(['Number of patients with at least one diagnosis below 45: ' num2str(patientsBelow45)]);
disp(['Number of patients with at least one diagnosis between 45 and 64: ' num2str(patients45to64)]);
disp(['Number of patients with at least one diagnosis at 65 or above: ' num2str(patients65AndAbove)]);
disp(['Percentage of patients with diagnosis below 45: ' num2str(patientsBelow45/totalPatients*100) '%']);
disp(['Percentage of patients with diagnosis between 45 and 64: ' num2str(patients45to64/totalPatients*100) '%']);
disp(['Percentage of patients with diagnosis at 65 or above: ' num2str(patients65AndAbove/totalPatients*100) '%']);

% Perform odds ratio analysis for all groups
resultsBelow45 = perform_odds_ratio_analysis(diagnosesBelow45, conditions, 'below45', analysisFolderPath);
results45to64 = perform_odds_ratio_analysis(diagnoses45to64, conditions, '45to64', analysisFolderPath);
results65AndAbove = perform_odds_ratio_analysis(diagnoses65AndAbove, conditions, '65plus', analysisFolderPath);

% Function to perform odds ratio analysis
function results = perform_odds_ratio_analysis(data, conditions, group_name, analysisFolderPath)
    % Note: 'data' contains only the diagnoses for the specific age group
    results = table('Size', [0, 19], ...
        'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'string', 'string', 'string', 'double', 'string', 'double', 'string', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'ConditionA', 'ConditionB', 'OddsRatio', 'PValue', 'OddsRatio_CI_Low', 'OddsRatio_CI_High', 'PairFrequency', 'Direction', 'Significance', 'Interpretation', 'AdjustedPValue', 'FDRSignificance', 'DirectionalFrequency', 'Precedence', 'DirectionalPercentage', 'Percentage', 'MedianDurationDays', 'MedianDurationYears', 'MedianDurationYearsWithIQR'});

    % Calculate the total number of patients with at least one diagnosis in this age group
    diagnosisCols = data.Properties.VariableNames(1:end-1);  % Exclude the YearOfBirth column
    totalPatientsInGroup = sum(any(~cellfun(@isempty, table2cell(data(:, diagnosisCols))), 2));

    for i = 1:length(conditions)
        for j = (i+1):length(conditions)
            cond1 = conditions{i};
            cond2 = conditions{j};
            
            [or, p, ci, pair_freq, direction, precedence, dir_freq, dir_percentage, median_duration_days, median_duration_years, median_duration_years_with_iqr] = compute_odds_ratio(data, cond1, cond2);
            
            significance = 'Not Significant';
            if p < 0.05
                significance = 'Significant';
            end
            
            percentage = (pair_freq / totalPatientsInGroup) * 100;

            newRow = table({cond1}, {cond2}, or, p, ci(1), ci(2), pair_freq, {direction}, {significance}, {''}, NaN, {''}, dir_freq, {precedence}, dir_percentage, percentage, median_duration_days, median_duration_years, {median_duration_years_with_iqr}, ...
                'VariableNames', {'ConditionA', 'ConditionB', 'OddsRatio', 'PValue', 'OddsRatio_CI_Low', 'OddsRatio_CI_High', 'PairFrequency', 'Direction', 'Significance', 'Interpretation', 'AdjustedPValue', 'FDRSignificance', 'DirectionalFrequency', 'Precedence', 'DirectionalPercentage', 'Percentage', 'MedianDurationDays', 'MedianDurationYears', 'MedianDurationYearsWithIQR'});
            results = [results; newRow];
        end
    end

    % Adjust p-values for multiple comparisons (FDR)
    [adjusted_pvals, fdr_significance] = fdr_bh(results.PValue);
    results.AdjustedPValue = adjusted_pvals;
    results.FDRSignificance = cellstr(string(fdr_significance));

    % Interpretation
    results.Interpretation = cell(height(results), 1);
    for i = 1:height(results)
        if results.AdjustedPValue(i) < 0.05
            results.Interpretation{i} = 'Significant after FDR correction';
        else
            results.Interpretation{i} = 'Not significant after FDR correction';
        end
    end

    % Save all results
    save(fullfile(analysisFolderPath, ['directional_odds_ratio_analysis_' group_name '.mat']), 'results', 'totalPatientsInGroup');
    
    % Include totalPatientsInGroup in the CSV file
    results.TotalPatientsInGroup = repmat(totalPatientsInGroup, height(results), 1);
    writetable(results, fullfile(analysisFolderPath, ['directional_odds_ratio_analysis_' group_name '.csv']));

    % Filter and save only FDR-significant associations with OR >= 2, CI_Low >= 2, and PairFrequency >= 10
    fdr_significant_high_freq_results = results(results.AdjustedPValue < 0.05 & ...
                                                results.OddsRatio >= 2 & ...
                                                results.OddsRatio_CI_Low >= 2 & ...
                                                results.PairFrequency >= 70, :);
    
    % Sort the results by odds ratio in descending order
    fdr_significant_high_freq_results = sortrows(fdr_significant_high_freq_results, 'OddsRatio', 'descend');
    
    save(fullfile(analysisFolderPath, ['fdr_significant_high_freq_odds_ratio_analysis_' group_name '.mat']), 'fdr_significant_high_freq_results', 'totalPatientsInGroup');
    writetable(fdr_significant_high_freq_results, fullfile(analysisFolderPath, ['fdr_significant_high_freq_odds_ratio_analysis_' group_name '.csv']));

    disp(['Analysis completed and results saved for ' group_name ' group.']);
    disp(['Total patients with at least one diagnosis in this group: ' num2str(totalPatientsInGroup)]);
end

% Function to compute odds ratio for a pair of conditions
function [or, p, ci, pair_freq, direction, precedence, dir_freq, dir_percentage, median_duration_days, median_duration_years, median_duration_years_with_iqr] = compute_odds_ratio(data, cond1, cond2)
    % Convert date strings to datetime objects
    dates1 = datetime(data.(cond1), 'InputFormat', 'dd/MM/yyyy');
    dates2 = datetime(data.(cond2), 'InputFormat', 'dd/MM/yyyy');
    
    % Create binary vectors
    binary1 = ~isnat(dates1);
    binary2 = ~isnat(dates2);
    
    % Create contingency table
    a = sum(binary1 & binary2);
    b = sum(binary1 & ~binary2);
    c = sum(~binary1 & binary2);
    d = sum(~binary1 & ~binary2);
    
    % Calculate odds ratio and confidence interval
    [~, p, stats] = fishertest([a b; c d]);
    or = stats.OddsRatio;
    ci = stats.ConfidenceInterval;
    
    % Calculate pair frequency
    pair_freq = a;
    
    % Determine directionality
    if or > 1
        direction = 'Positive';
    elseif or < 1
        direction = 'Negative';
    else
        direction = 'Neutral';
    end
    
    % Calculate precedence and directional frequency
    earlier = sum((binary1 & binary2) & (dates1 < dates2));
    later = sum((binary1 & binary2) & (dates1 > dates2));
    
    if earlier > later
        precedence = [cond1 ' precedes ' cond2];
        dir_freq = earlier;
    elseif later > earlier
        precedence = [cond2 ' precedes ' cond1];
        dir_freq = later;
    else
        precedence = 'No clear precedence';
        dir_freq = a;
    end
    
    % Calculate directional percentage
    dir_percentage = (dir_freq / pair_freq) * 100;
    
    % Calculate median duration and IQR
    [median_duration_days, median_duration_years, median_duration_years_with_iqr] = calculate_median_duration_with_iqr(dates1, dates2);
end

% Simplified FDR correction function
function [adjusted_pvals, significance] = fdr_bh(pvals, q)
    if nargin < 2
        q = 0.05;
    end
    
    % Sort p-values
    [sorted_pvals, sortIdx] = sort(pvals);
    n = length(pvals);
    
    % Calculate FDR threshold
    thresholds = (1:n) * q / n;
    below_threshold = sorted_pvals <= thresholds';
    
    % Find the largest p-value that is below the threshold
    max_below = find(below_threshold, 1, 'last');
    
    if isempty(max_below)
        adjusted_pvals = ones(size(pvals));
    else
        adjusted_pvals = min(1, pvals * n / (1:n)');
    end
    

% Adjusted p-values
    adjusted_pvals = min(1, adjusted_pvals(sortIdx));
    
    % Significance based on FDR
    significance = adjusted_pvals < q;
end

% Function to calculate median duration between conditions with IQR
function [median_duration_days, median_duration_years, median_duration_years_with_iqr] = calculate_median_duration_with_iqr(dates1, dates2)
    valid_pairs = ~isnat(dates1) & ~isnat(dates2);
    durations = days(abs(dates2(valid_pairs) - dates1(valid_pairs)));
    median_duration_days = median(durations);
    median_duration_years = median_duration_days / 365.25; % Approximate conversion to years
    
    % Calculate IQR
    q1 = prctile(durations, 25) / 365.25;
    q3 = prctile(durations, 75) / 365.25;
    
    % Format the string with median and IQR
    median_duration_years_with_iqr = sprintf('%.2f [%.2f-%.2f]', median_duration_years, q1, q3);
end