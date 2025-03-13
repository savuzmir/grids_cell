function [out1, out2, out3] = computePsychometricCurve(left_mag, left_prob, right_mag, right_prob, chosen)
% takes in structure with task data (inp) that can be split by brainer 
% returns probability of choice, and EV between left and right

% chosen is 1 for right  and vice versa

% discretized
leftEV  = discretize(left_mag, 5).*discretize(left_prob, 5);
rightEV = discretize(right_mag, 5).*discretize(right_prob, 5);

EV_DIFF = rightEV-leftEV;

unqEV = unique(EV_DIFF);

out3 =  unqEV;

out1 = [];
out2 = [];

for i = 1:length(unqEV) % we define a static number for now for the info gathering dataset 

    relTr = find(EV_DIFF==unqEV(i));
    
    out1(i, :) = sum(chosen(relTr))/length(relTr);
    out2(i, :) = length(relTr);
end



end

