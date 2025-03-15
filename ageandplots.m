clear
% Specify the folder paths
baseFolder = 'C:\Users\cogc\OneDrive - Loughborough University\MATLAB\MATLAB\ALLDataDatesHere\PERFECTCODE\SplitInYear';
maleAnalysisFolderPath = fullfile(baseFolder, 'males');
femaleAnalysisFolderPath = fullfile(baseFolder, 'females');
dataFolder = baseFolder;

% Create the analysis folders if they don't exist
for folderPath = {maleAnalysisFolderPath, femaleAnalysisFolderPath}
    if ~exist(folderPath{1}, 'dir')
        mkdir(folderPath{1});
        disp(['Created analysis folder: ', folderPath{1}]);
    else
        disp(['Analysis folder already exists: ', folderPath{1}]);
    end
end

% Load the data
load(fullfile(dataFolder, 'maledatebirth.mat'), 'datamerge');
maleDateMerge = datamerge;
load(fullfile(dataFolder, 'femaledatebirth.mat'), 'datamerge');
femaleDateMerge = datamerge;

% Get all unique condition names from both datasets
maleConditions = maleDateMerge.Properties.VariableNames(1:end-1);
femaleConditions = femaleDateMerge.Properties.VariableNames(1:end-1);
allConditions = unique([maleConditions, femaleConditions]);

% Identify conditions only present in female dataset
femaleOnlyConditions = setdiff(femaleConditions, maleConditions);

% Create a set of conditions for male plot
maleRelevantConditions = setdiff(allConditions, femaleOnlyConditions);

% Function to process data and calculate ages
function ages_table = processData(datamerge, allConditions)
    ages_at_diagnosis = cell(height(datamerge), length(allConditions));
    dataConditions = datamerge.Properties.VariableNames(1:end-1);

    % Process each patient
    for i = 1:height(datamerge)
        yearOfBirth = datamerge{i, end};
        
        for j = 1:length(allConditions)
            if ismember(allConditions{j}, dataConditions)
                diagnosisDate = datamerge{i, allConditions{j}}{1};  % Access the string inside the cell
                if ~isempty(diagnosisDate) && ~isempty(strtrim(diagnosisDate))
                    diagDate = datetime(diagnosisDate, 'InputFormat', 'dd/MM/yyyy');
                    ageAtDiagnosis = year(diagDate) - yearOfBirth;
                    ages_at_diagnosis{i, j} = ageAtDiagnosis;
                else
                    ages_at_diagnosis{i, j} = NaN;
                end
            else
                ages_at_diagnosis{i, j} = NaN;
            end
        end
    end

    % Convert cell array to table
    ages_table = cell2table(ages_at_diagnosis, 'VariableNames', allConditions);
    ages_table.YearOfBirth = datamerge{:, end};
end

% Process male and female data
maleAgesTable = processData(maleDateMerge, allConditions);
femaleAgesTable = processData(femaleDateMerge, allConditions);

% Save results
save(fullfile(maleAnalysisFolderPath, 'ages_at_diagnosis.mat'), 'maleAgesTable');
writetable(maleAgesTable, fullfile(maleAnalysisFolderPath, 'ages_at_diagnosis.csv'));
save(fullfile(femaleAnalysisFolderPath, 'ages_at_diagnosis.mat'), 'femaleAgesTable');
writetable(femaleAgesTable, fullfile(femaleAnalysisFolderPath, 'ages_at_diagnosis.csv'));

% Prepare data for box plot
maleAgesForPlot = table2array(maleAgesTable(:, 1:end-1));
femaleAgesForPlot = table2array(femaleAgesTable(:, 1:end-1));

% Calculate median ages
femaleMedianAges = nanmedian(femaleAgesForPlot);
maleMedianAges = nanmedian(maleAgesForPlot);

greenColor = [0, 0.5, 0]; % Dark green
% Create separate box plots for females (top) and males (bottom)
figure('Position', [100, 100, 1500, 1000]);

% Female box plot (top)
subplot(2,1,1);
boxplot(femaleAgesForPlot, 'Labels', allConditions, 'Colors', greenColor);
hold on;
plot(1:length(allConditions), femaleMedianAges, 'r*', 'MarkerSize', 5);
for i = 1:length(allConditions)
    text(i, femaleMedianAges(i), sprintf('%.1f', femaleMedianAges(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'black');
end
hold off;
title({'Female Age at Diagnosis by Condition', '\newline'});
ylabel('Age');
ylim([0, 100]); % Adjust y-axis limits as needed
xtickangle(45); % Rotate x-axis labels for better readability
set(gca, 'TickLength', [0 0]); % Remove tick lines

% Male box plot (bottom)
subplot(2,1,2);
maleRelevantAgesForPlot = maleAgesForPlot(:, ismember(allConditions, maleRelevantConditions));
maleRelevantMedianAges = maleMedianAges(ismember(allConditions, maleRelevantConditions));
boxplot(maleRelevantAgesForPlot, 'Labels', maleRelevantConditions, 'Colors', greenColor);
hold on;
plot(1:length(maleRelevantConditions), maleRelevantMedianAges, 'r*', 'MarkerSize', 5);
for i = 1:length(maleRelevantConditions)
    text(i, maleRelevantMedianAges(i), sprintf('%.1f', maleRelevantMedianAges(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'black');
end
hold off;
title('Male Age at Diagnosis by Condition');
ylabel('Age');
ylim([0, 100]);  % Adjust y-axis limits as needed
xtickangle(45);  % Rotate x-axis labels for better readability
set(gca, 'TickLength', [0 0]); % Remove tick lines

% Save the plot
saveas(gcf, fullfile(baseFolder, 'age_at_diagnosis_boxplot_by_gender.png'));
saveas(gcf, fullfile(baseFolder, 'age_at_diagnosis_boxplot_by_gender.fig'));

% Calculate summary statistics
calculateSummaryStats = @(agesForPlot) ...
    [nanmean(agesForPlot); nanmedian(agesForPlot); nanstd(agesForPlot); ...
     sum(~isnan(agesForPlot)); ...
     sum(~isnan(agesForPlot)) ./ size(agesForPlot, 1) * 100];

maleSummary = calculateSummaryStats(maleAgesForPlot);
femaleSummary = calculateSummaryStats(femaleAgesForPlot);

% Calculate overall statistics
maleOverallMedian = nanmedian(maleAgesForPlot(:));
femaleOverallMedian = nanmedian(femaleAgesForPlot(:));
maleMedianRange = [min(maleSummary(2,:)), max(maleSummary(2,:))];
femaleMedianRange = [min(femaleSummary(2,:)), max(femaleSummary(2,:))];

% Combine male and female statistics
combinedStats = array2table([maleSummary; femaleSummary], ...
    'VariableNames', allConditions, ...
    'RowNames', {'Male_Mean', 'Male_Median', 'Male_StdDev', 'Male_Count', 'Male_Percentage', ...
                 'Female_Mean', 'Female_Median', 'Female_StdDev', 'Female_Count', 'Female_Percentage'});

% Add overall statistics
overallStats = table(maleOverallMedian, femaleOverallMedian, maleMedianRange(1), maleMedianRange(2), femaleMedianRange(1), femaleMedianRange(2), ...
    'VariableNames', {'Male_Overall_Median', 'Female_Overall_Median', 'Male_Median_Min', 'Male_Median_Max', 'Female_Median_Min', 'Female_Median_Max'}, ...
    'RowNames', {'Value'});

% Save combined statistics
writetable(combinedStats, fullfile(baseFolder, 'combined_age_statistics.csv'), 'WriteRowNames', true);
writetable(overallStats, fullfile(baseFolder, 'overall_age_statistics.csv'), 'WriteRowNames', true);

disp('Combined age statistics saved to combined_age_statistics.csv');
disp('Overall age statistics saved to overall_age_statistics.csv');

% Display overall statistics
disp('Overall Statistics:');
disp(overallStats);

% Display summary for males
totalMalePatients = height(maleAgesTable);
malePatientsWithDiagnosis = sum(any(~isnan(maleAgesForPlot), 2));
disp(['Total number of male patients: ' num2str(totalMalePatients)]);
disp(['Number of male patients with at least one diagnosis: ' num2str(malePatientsWithDiagnosis)]);
disp(['Percentage of male patients with diagnosis: ' num2str(malePatientsWithDiagnosis/totalMalePatients*100) '%']);

% Display summary for females
totalFemalePatients = height(femaleAgesTable);
femalePatientsWithDiagnosis = sum(any(~isnan(femaleAgesForPlot), 2));
disp(['Total number of female patients: ' num2str(totalFemalePatients)]);
disp(['Number of female patients with at least one diagnosis: ' num2str(femalePatientsWithDiagnosis)]);
disp(['Percentage of female patients with diagnosis: ' num2str(femalePatientsWithDiagnosis/totalFemalePatients*100) '%']);