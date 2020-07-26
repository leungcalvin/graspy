function fname = separateED(directory,stateID,DHForCI,varargin)
%%% Input
% varargin{1} = suffix
%%% output
% fname: names of files generated.

% stateID = '2p_3';
% directory = './';
% DHForCI = 'CI';

fname = {};

if strcmp(DHForCI,'DHF')
    fname_load = [stateID,'.ed'];
elseif strcmp(DHForCI,'CI')
    fname_load = [stateID,'.ced'];
else
    error('Invalid input of DHForCI.')
end

%% load every line in file
fstr = {};
i = 0;
fid = fopen([directory,fname_load]);
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

%% save each level in its own file
% read NNNP
str_NNNP = fstr{2};
expr_NNNP = '\sNNNP =\s+(?<NNNP>\d+)';
% expr_csf = ['(?<subshell>\d+',expr_lterm,'-?)\s?\(\s?(?<occupn>\d+)\)'];
S_NNNP = regexp(str_NNNP,expr_NNNP,'names');
NNNP = str2double(S_NNNP.NNNP);

% number of level
N = (n_line-2-1)/(NNNP+4);

% suffix
if length(varargin) >= 1
    if length(varargin{1}) ~= N
    error('Length of suffix is not the same to the number of levels.')
    end
    fname_suffix = varargin{1};
else
    fname_suffix = compose(['_No%0',num2str(floor(N/10)+1),'i'],1:N);
end

for No = 1:N
%     fname_save = [directory,stateID,'_level',...
%         sprintf(['%0',num2str(floor(N_level/10)+1),'i'],level)];
    fname_save = [stateID,fname_suffix{No}];
    if strcmp(DHForCI,'DHF')
        fname_save = [fname_save,'.ed'];
    else
        fname_save = [fname_save,'.ced'];
    end
    
    fname{No} = fname_save;
    fid = fopen([directory,fname_save],'w');
    
    str_start = fstr(1:2);
    str_start{1} = sprintf('%s, level %i in %s',...
        str_start{1},No,fname_load);
    fprintf(fid,'%s\n',fstr{1:2});
    for li = (2+(NNNP+4)*(No-1)+1):(2+(NNNP+4)*No)
        fprintf(fid,'%s\n',fstr{li});
    end
    
    str_end = fstr{end};
    fprintf(fid,str_end);
    
    fclose(fid);
end

fprintf('files for level 1 - %i generated.\n',N)
end