function sessIndx = extract_sess(nSessions, trialsSliced)

%       we want to find IDs of the first occurence of a given session to
%       extract the corresponding data
        sessIndx = [];
        for i = 1:size(nSessions, 1)
            tmp = trialsSliced(nSessions(i) == trialsSliced(:, 2), 1);
            sessIndx(i) = tmp(1);
        end
        
end