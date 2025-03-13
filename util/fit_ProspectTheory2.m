function [energy, optim_param] = fit_ProspectTheory2(AuxiliaryLFP, sessIndx, sess, startingParam)

parameters = startingParam;

    function [energy] = fit_ProspectTheorySingle(parameters)
        
        alpha = parameters(1);
        gamma = parameters(2);
        theta = parameters(3);
        delta = parameters(4);
        
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
      
        % get choice
        choice = (AuxiliaryLFP(sessIndx(sess)).Right_pay == AuxiliaryLFP(sessIndx(sess)).ChMR) + 1;

        % noise
        probL = (1-delta).*(1./(1+exp(-theta*(valLeft - valRight)))) + delta*0.5; 
        probR = (1-delta).*(1./(1+exp(-theta*(valRight - valLeft)))) + delta*0.5; 

        
        probPT = [probL, probR];
        
        energy = [];
        for tr = 1:size(choice, 1)
            energy(tr) = [probPT(tr, choice(tr))];
        end
        energy = -sum(log(energy));
    end


lb = [10e-2, 10e-2, 10e-2, 10e-2];
ub = [1, 1, 10, 1];

A = [];
b = [];
Aeq = [];
beq = [];

[energy, optim_param] = fmincon(@fit_ProspectTheorySingle, parameters, A, b, Aeq, beq, lb, ub, []);

end