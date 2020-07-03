function T_level = readRlevels(fname)
% fname = './rlevelseV.out';

%% load every line in file
fstr = {};
i = 0;
fid = fopen(fname);
while ~feof(fid)
    i = i+1;
    fstr{i} = fgetl(fid);
    %     disp(fstr{i})
end
fclose(fid);
fstr = fstr.';

N_line = length(fstr);

% % check if it is right file
% if ~strcmp(fstr{1},'G92RWF')
%     error('wrong .w.readrwf file.')
% end

%% read data
% find starting line
headerstr = ' No Pos  J';
for li = 1:N_line
    if ~isempty(strfind(fstr{li},headerstr))
        li_header = li;
        break;
    end
end
li_start = li_header + 3;

endstr = '----';
for li = li_start:N_line
    if ~isempty(strfind(fstr{li},endstr))
        li_end = li - 1;
        break;
    end
end

% read out numbers
expr = '^\s*\d+\s+\d+\s+(?<J>\d+(/\d+)?)\s+(?<Parity>(+|-)?)\s+(?<E_tot>-?\d+\.\d+)\s+(?<Level>-?\d+\.\d+)\s+(?<Spliting>-?\d+\.\d+)(\s+(?<Configuration>.+))?';
% C = {};
% lvi = 0;
% for li = li_start:N_line
%     lvi = lvi + 1;
%     C{lvi} = regexp(fstr{li},expr,'tokens');
%     if isempty(C{lvi})
%         C = C(1:end-1);
%         break;
%     end
% end
%
% C = vertcat(C{:});
% C = vertcat(C{:});

N_level = li_end - li_start + 1;
lvi = 0;
for li = li_start:li_end
    lvi = lvi + 1;
    S = regexp(fstr{li},expr,'names');
    if isempty(S)
        continue;
    end
    S_level(lvi) = S;
end
T_level = struct2table(S_level);
T_level = [T_level(:,end),T_level(:,1:end-1)];

% T_level.J = str2double(T_level.J);
T_level.E_tot = str2double(T_level.E_tot);
T_level.Level = str2double(T_level.Level);
T_level.Spliting = str2double(T_level.Spliting);

end