    function [energy, optim_param] = fit_ProspectTheory3(AuxiliaryLFP, sessIndx, sess, startingParam)

parameters = startingParam;

    function [energy] = fit_ProspectTheorySingle(parameters)
        
        alpha = parameters(1);
        gamma = parameters(2);
        theta = parameters(3);
        delta = parameters(4);
        zeta  = parameters(5);
        
        nTrials = size(AuxiliaryLFP(sessIndx(sess)).Prob_top, 2);
        
        rawMagL = AuxiliaryLFP(sessIndx(sess)).Left_pay;
        rawMagR = AuxiliaryLFP(sessIndx(sess)).Right_pay;
        
        rawProbL = AuxiliaryLFP(sessIndx(sess)).Left_prob;
        rawProbR = AuxiliaryLFP(sessIndx(sess)).Right_prob;
        
        magL = rawMagL.^alpha;                         
        probL = exp(-(-log(rawProbL)).^gamma);
        
        magR = rawMagR.^alpha;
        probR = exp(-(-log(rawProbR)).^gamma);
        
        valLeft = magL.*probL;                          
        valRight = magR.*probR;

        % chosen option - right
        choice = (AuxiliaryLFP(sessIndx(sess)).Right_pay == AuxiliaryLFP(sessIndx(sess)).ChMR) + 1;
        
        % 1 = same, 0 = not same
        perseveranceL = [1; abs(diff(choice))]; 
        perseveranceR = [1; abs(diff(choice))]; 
      
        % correct for 1st
        perseveranceL(1) = 0;
        perseveranceR(1) = 0;
        
        for i = 2:nTrials
            perseveranceL(i) = (choice(i-1) == choice(i)) & (choice(i) == 1);
            perseveranceR(i) = (choice(i-1) == choice(i)) & (choice(i) == 2);
        end        

        % % adding noise
        probL = (1-delta).*(1./(1+exp(-theta*(valLeft - valRight + perseveranceL*zeta)))) + delta*(1/(1+exp(perseveranceL*zeta)))'; 
        probR = (1-delta).*(1./(1+exp(-theta*(valRight - valLeft + perseveranceR*zeta)))) + delta*(1/(1+exp(perseveranceR*zeta)))'; 

        probPT = [probL, probR];
        
        energy = [];
        for tr = 1:size(choice, 1)
            energy(tr) = [probPT(tr, choice(tr))];
        end
        energy = -sum(log(energy));
    end


lb = [10e-2, 10e-2, 10e-2, 10e-2, -10e-2];
ub = [1, 1, 10, 1, 10e-2];

A = [];
b = [];
Aeq = [];
beq = [];

[energy, optim_param] = fmincon(@fit_ProspectTheorySingle, parameters, A, b, Aeq, beq, lb, ub, []);

end