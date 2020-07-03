%%% Read .(c)ed file from REDF routine in GRASP2018
function [r,ED] = readED(directory,stateID,DHForCI)

% stateID = '2s_3';
% directory = './';
% DHForCI = 'CI';

if strcmp(DHForCI,'DHF')
    fname_load = [directory,stateID,'.ed'];
elseif strcmp(DHForCI,'CI')
    fname_load = [directory,stateID,'.ced'];
else
    error('Invalid input of DHForCI.')
end

%% load every line in file
fstr = {};
i = 0;
fid = fopen(fname_load);
while ~feof(fid)
    i = i+1;
    fstr{i} = fgetl(fid);
    %     disp(fstr{i})
end
fclose(fid);
fstr = fstr.';

n_line = length(fstr);

% check if it is right file
str_chk = ' REDF1: Radial (volume) electron density function over r';
if ~contains(fstr{1},str_chk)
    error('wrong .(c)ed file.')
end
str_chk = ' REDF1: End of file';
if ~strcmp(fstr{end},str_chk)
    error('wrong .(c)ed file.')
end

%% Read data
% read NNNP
str_NNNP = fstr{2};
expr_NNNP = '^\sNNNP =\s+(?<NNNP>\d+)';
% expr_csf = ['(?<subshell>\d+',expr_lterm,'-?)\s?\(\s?(?<occupn>\d+)\)'];
S_NNNP = regexp(str_NNNP,expr_NNNP,'names');
NNNP = str2double(S_NNNP.NNNP);

% read r, ED
datastr = fstr(6:NNNP+6-1);
datastr = cellfun(@(c) replace(c,'D','E'),datastr,'UniformOutput',false);
data = cellfun(@(c) [str2double(c(1:23)),str2double(c(24:end))],datastr,...
    'UniformOutput',false);
data = vertcat(data{:});
r = data(:,1); % a_0; position r
ED = data(:,2); % a_0^{-3}; (volume) radial electron number density

% % test plot
% test_fig = figure;
% hold on;
% ax = test_fig.CurrentAxes;
% h = plot(r,ED);
% ax.XScale = 'log';
% ax.XLim = r([1,end]);
% xlabel('r (a_0)')
% ylabel('Number Density (a_0^{-3})')

end
