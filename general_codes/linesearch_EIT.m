function [step,F,Uref] = linesearch_EIT(x,d,F0,Ffunc,firststep,minval,nmax,figno,DRAW,varargin);
%
% x = previous estimate
% d = search direction
% F0 = functional value corresp. to x
% Ffunc = function used for computing the functional to be minimized (string)

%keyboard


% FIRST STEP



stepvec = [0 firststep];
xnew = x + firststep*d;
negind = find(xnew < minval);
xnew(negind) = minval;
Fstr = [Ffunc '(xnew'];
for ii=1:length(varargin)
  eval(sprintf('P%d = varargin{%d};',ii,ii));
  Fstr = [Fstr sprintf(',P%d',ii)];
end
Fstr = [Fstr ');'];
%disp(['F-Call: ' Fstr])

[Fnew,Urefnew] = eval(Fstr); % compute the functional at step(iter)
F_ALL = [F0 Fnew];
%Uref_ALL = [Uref0,Urefnew]; 
[Fmin,Fminind] = min(F_ALL);
n = 1;

if DRAW
  figure(figno), clf, 
  plot(stepvec,F_ALL,'b+')
  title('Line-search'), drawnow
end

% FIND THE LIMIT

if(Fminind == 2)   % DOUBLE THE STEP LENGTH REPETITIVELY

  while(Fminind == length(F_ALL) && n < nmax)
    newstep = 2*stepvec(end);
    stepvec = [stepvec newstep];
    xnew = x + newstep*d;
    negind = find(xnew < minval);
    xnew(negind) = minval;
    Fstr = [Ffunc '(xnew'];
    for ii=1:length(varargin)
       eval(sprintf('P%d = varargin{%d};',ii,ii));
       Fstr = [Fstr sprintf(',P%d',ii)];
    end
    Fstr = [Fstr ');'];
    %disp(['F-Call: ' Fstr])
    [Fnew,Urefnew] = eval(Fstr); % compute the functional at step(iter)
    F_ALL = [F_ALL Fnew];
    [Fmin,Fminind] = min(F_ALL);
    %Uref_ALL = [Uref_ALL, Urefnew];
    n = n+1;
    if DRAW
      figure(figno), clf, 
      plot(stepvec,F_ALL,'b+')
      hold on
      plot(newstep,Fnew,'ro')
      title('Line-search'), drawnow
    end
  end

else  % HALF THE STEP LENGTH REPETITIVELY

  while(Fminind == 1 && n<nmax)
    newstep = .5*stepvec(2);
    stepvec = [stepvec(1) newstep stepvec(2:end)];
    xnew = x + newstep*d;
    negind = find(xnew < minval);
    xnew(negind) = minval;
    Fstr = [Ffunc '(xnew'];
    for ii=1:length(varargin)
       eval(sprintf('P%d = varargin{%d};',ii,ii));
       Fstr = [Fstr sprintf(',P%d',ii)];
    end
    Fstr = [Fstr ');'];
    %disp(['F-Call: ' Fstr])
    [Fnew,Urefnew] = eval(Fstr); % compute the functional at step(iter)
    F_ALL = [F_ALL(1) Fnew F_ALL(2:end)];
    [Fmin,Fminind] = min(F_ALL);
    %Uref_ALL = [Uref_ALL(1) Urefnew Uref_ALL(2:end)];
    n = n+1;
    if DRAW
      figure(figno), clf, 
      plot(stepvec,F_ALL,'b+')
      hold on
      plot(newstep,Fnew,'ro')
      title('Line-search'), drawnow
    end
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
    xnew = x + newstep*d;
    negind = find(xnew < minval);
    xnew(negind) = minval;
    Fstr = [Ffunc '(xnew'];
    for ii=1:length(varargin)
       eval(sprintf('P%d = varargin{%d};',ii,ii));
       Fstr = [Fstr sprintf(',P%d',ii)];
    end
    Fstr = [Fstr ');'];
    %disp(['F-Call: ' Fstr])
    [Fnew,Urefnew] = eval(Fstr); % compute the functional at step(iter)
    F_ALL = [F_ALL(1:Fminind) Fnew F_ALL(Fminind+1:end)];
    %Uref_ALL = [Uref_ALL(1:Fminind) Urefnew Uref_ALL(Fminind+1:end)];
    [Fmin,Fminind] = min(F_ALL);
    if DRAW
      figure(figno), clf, 
      plot(stepvec,F_ALL,'b+')
      hold on
      plot(newstep,Fnew,'ro')
      title('Line-search'), drawnow
    end
    n = n+1;

    if (Fnew > Fmin)

      % SMALLER STEP

      newstep = mean(stepvec([Fminind-1 Fminind]));
      stepvec = [stepvec(1:Fminind-1) newstep stepvec(Fminind:end)];
      xnew = x + newstep*d;
      negind = find(xnew < minval);
      xnew(negind) = minval;
      Fstr = [Ffunc '(xnew'];
      for ii=1:length(varargin)
         eval(sprintf('P%d = varargin{%d};',ii,ii));
         Fstr = [Fstr sprintf(',P%d',ii)];
      end
      Fstr = [Fstr ');'];
      %disp(['F-Call: ' Fstr])
      [Fnew,Urefnew] = eval(Fstr); % compute the functional at step(iter)
      F_ALL = [F_ALL(1:Fminind-1) Fnew F_ALL(Fminind:end)];
      %Uref_ALL = [Uref_ALL(1:Fminind-1) Urefnew Uref_ALL(Fminind:end)];
      [Fmin,Fminind] = min(F_ALL);
      if DRAW
        figure(figno), clf, 
        plot(stepvec,F_ALL,'b+')
        hold on
        plot(newstep,Fnew,'ro')
        title('Line-search'), drawnow
      end
      n = n+1;
    end

  end

  step = stepvec(Fminind);
  F = F_ALL(Fminind);
  %Uref = Uref_ALL(Fminind);
  Uref = [];

end

