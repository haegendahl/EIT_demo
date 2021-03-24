function [step,F,Uref]  = linesearch2(FIRSTSTEP,nmax, DIRECTION, Uel, F0, Node, Element, I, THETA, P1st, z, MeasPatt, Ln, R, Ai, beta, alpha, figno)

% change to 1st order fwdsolution


stepvec = [0 FIRSTSTEP];

thetanew = THETA + FIRSTSTEP*DIRECTION;

%tmp = sigmaest+steps(i)*stepdirection;

negind = thetanew <= 0;
thetanew(negind) = 10^(-6);


% Evaluate the Functional || x || at the FIRST iteration
% compute the functional at step(iter)

adm = thetanew;

Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*adm,z,MeasPatt,'real');
Urefel = Urefel(:);

Fnew = 0.5*norm(Ln*(Uel-Urefel))^2 +alpha*TV3D_functional(R, Ai, thetanew, beta);



%----------------------END FIRST ITERATION ----------------------


F_ALL = [F0 Fnew];

[Fmin, Fminind] = min(F_ALL);
n = 1;


figure(figno), clf
plot(stepvec,F_ALL,'b+')
drawnow


% DOUBLE THE STEP LENGTH REPETITIVELY

if(Fminind == 2)
    
    while(Fminind == length(F_ALL) && n < nmax)
        newstep = 2*stepvec(end);
        stepvec = [stepvec newstep];
        
        thetanew = THETA + newstep*DIRECTION;
        %xnew = x + newstep*d;
        negind = thetanew <= 0;
        thetanew(negind) = 10^(-6);
        
        
        
        %[Fnew,Urefnew] = eval(Fstr);  compute the functional at step(iter)
        % compute the functional at step(iter)
        
     
        adm = thetanew;
        
        Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*adm,z,MeasPatt,'real');
        Urefel = Urefel(:);
        
        Fnew = 0.5*norm(Ln*(Uel-Urefel))^2 + alpha*TV3D_functional(R, Ai, thetanew, beta);
        
        F_ALL = [F_ALL Fnew];
        
        [Fmin,Fminind] = min(F_ALL);
        %Uref_ALL = [Uref_ALL, Urefnew];
        n = n+1;
        figure(figno), clf
        plot(stepvec,F_ALL,'b+')
        hold on
        plot(newstep,Fnew,'ro')
        drawnow
        
    end
    
else  % HALF THE STEP LENGTH REPETITIVELY
    
    while(Fminind == 1 && n<nmax)
        newstep = .5*stepvec(2);
        stepvec = [stepvec(1) newstep stepvec(2:end)];
        thetanew = THETA + newstep*DIRECTION;
        negind = thetanew <= 0;
        thetanew(negind) = 10^(-6);
        
        %[Fnew,Urefnew] = eval(Fstr);  compute the functional at step(iter)
        % compute the functional at step(iter)
        
        adm = thetanew;
        
        Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*adm,z,MeasPatt,'real');
        Urefel = Urefel(:);
        
        Fnew = 0.5*norm(Ln*(Uel-Urefel))^2 + alpha*TV3D_functional(R, Ai, thetanew, beta);
        
        F_ALL = [F_ALL(1) Fnew F_ALL(2:end)];
        
        [Fmin,Fminind] = min(F_ALL);
        %Uref_ALL = [Uref_ALL(1) Urefnew Uref_ALL(2:end)];
        n = n+1;
        figure(figno), clf
        
        plot(stepvec,F_ALL,'b+')
        hold on
        plot(newstep,Fnew,'ro')
        drawnow
    end
    
end


if (Fminind == 1) % GIVE UP
    
    step = 0;
    F = F0;
    %Uref = Uref0;
    Uref = [];
    
else
    
    while(n < nmax)
        
        % LARGER STEP
        
        newstep = mean(stepvec([Fminind Fminind+1]));
        stepvec = [stepvec(1:Fminind) newstep stepvec(Fminind+1:end)];
        thetanew = THETA + newstep*DIRECTION;
        
        negind = thetanew <= 0;
        thetanew(negind) = 10^(-6);
        
        % compute the functional at step(iter)
               
        adm = thetanew;
        
        Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*adm,z,MeasPatt,'real');
        Urefel = Urefel(:);
        
        Fnew = 0.5*norm(Ln*(Uel-Urefel))^2 + alpha*TV3D_functional(R, Ai, thetanew, beta);
        
        F_ALL = [F_ALL(1:Fminind) Fnew F_ALL(Fminind+1:end)];
        %Uref_ALL = [Uref_ALL(1:Fminind) Urefnew Uref_ALL(Fminind+1:end)];
        
        [Fmin,Fminind] = min(F_ALL);
        
        figure(figno), clf
        plot(stepvec,F_ALL,'b+')
        hold on
        plot(newstep,Fnew,'ro')
        drawnow
        n = n+1;
        
        if (Fnew > Fmin)
            
            % SMALLER STEP
            
            newstep = mean(stepvec([Fminind-1 Fminind]));
            stepvec = [stepvec(1:Fminind-1) newstep stepvec(Fminind:end)];
            thetanew = THETA + newstep*DIRECTION;
            negind = thetanew <= 0;
            thetanew(negind) = 10^(-6);
            
            
            
            % compute the functional at step(iter)
            adm = thetanew;
            
            Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*adm,z,MeasPatt,'real');
            Urefel = Urefel(:);
            
            Fnew = 0.5*norm(Ln*(Uel-Urefel))^2 + alpha*TV3D_functional(R, Ai, thetanew, beta);
            
            
            F_ALL = [F_ALL(1:Fminind-1) Fnew F_ALL(Fminind:end)];
            
            [Fmin,Fminind] = min(F_ALL);
            figure(figno), clf
            plot(stepvec,F_ALL,'b+')
            hold on
            plot(newstep,Fnew,'ro')
            drawnow
            n = n+1;
        end
        
    end
    
    
    step = stepvec(Fminind);
    F = F_ALL(Fminind);
    %Uref = Uref_ALL(Fminind);
    Uref = [];
    
end