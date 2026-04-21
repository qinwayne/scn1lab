%% repeated anova test if stages share the same mean
% input = fish (35) x stages (20)
function [output,AT] = ranova2(input, varargin)
T = array2table(input);
T.fish = [1:size(input,1)]';
% alpha1 = 'abcdefghijklmnopqrstuvwxyz';
% alpha_cell_array = sprintfc('%c',alpha1);
% T.Properties.VariableNames = alpha_cell_array(1:size(input,2));
T.genotypes = {'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'WT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT' 'MUT'}';
% writetable(T,'C:\Users\wqin2\OneDrive - The University of Melbourne\Research\SCN1\stats\GT states\gt.csv')
% create the within-subjects design
withinDesign = table([1:size(input,2)]','VariableNames',{'Stage'});
withinDesign.Stage = categorical(withinDesign.Stage);

% create the repeated measures model and do the anova
rm = fitrm(T,['input1-input' mat2str(size(input,2)) ' ~ genotypes'],'WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Stage');
% AT = ranova(rm) % ,'WithinModel','Stage' remove comma to see ranova's table

disp(anovaTable(AT, 'Measure (units)'));
if nargin == 1
    output = multcompare(rm,'genotypes','By','Stage','ComparisonType','bonferroni');
elseif varargin{1} == 2
    output = multcompare(rm,'Stage','By','genotypes','ComparisonType','bonferroni');
elseif varargin{1} == 3
    output = multcompare(rm,'Stage','ComparisonType','bonferroni');
elseif varargin{1} == 4
    output = multcompare(rm,'genotypes','ComparisonType','bonferroni');
elseif varargin{1} == 5
    output = multcompare(rm,'Stage','By','Stage:genotypes','ComparisonType','bonferroni');
end

end


% -------------------------------------------------------------------------
% Function to create a conventional ANOVA table from the overly-complicated
% and confusing ANOVA table created by the ranova function.
function [s] = anovaTable(AT, dvName)
c = table2cell(AT);
% remove erroneous entries in F and p columns
for i=1:size(c,1)
    if c{i,4} == 1
        c(i,4) = {''};
    end
    if c{i,5} == .5
        c(i,5) = {''};
    end
end
% use conventional labels in Effect column
effect = AT.Properties.RowNames;
for i=1:length(effect)
    tmp = effect{i};
    tmp = erase(tmp, '(Intercept):');
    tmp = strrep(tmp, 'Error', 'Participant');
    effect(i) = {tmp};
end
% determine the required width of the table
fieldWidth1 = max(cellfun('length', effect)); % width of Effect column
fieldWidth2 = 57; % width for df, SS, MS, F, and p columns
barDouble = repmat('=', 1, fieldWidth1 + fieldWidth2);
barSingle = repmat('-', 1, fieldWidth1 + fieldWidth2);
% re-organize the data
c = c(2:end,[2 1 3 4 5]);
c = [num2cell(repmat(fieldWidth1, size(c,1), 1)), effect(2:end), c]';
% create the ANOVA table
s = sprintf('ANOVA table for %s\n', dvName);
s = [s sprintf('%s\n', barDouble)];
s = [s sprintf('%-*s %4s %11s %14s %9s %9s\n', fieldWidth1, 'Effect', 'df', 'SS', 'MS', 'F', 'p')];
s = [s sprintf('%s\n', barSingle)];
s = [s sprintf('%-*s %4d %14.5f %14.5f %10.3f %10.4f\n', c{:})];
s = [s sprintf('%s\n', barDouble)];
end