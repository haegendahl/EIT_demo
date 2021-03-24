%% THIS IS NOT A FUNCTION, EXECUTE WITH CARE !! %%

   if bmax>10, bmax=.5*bmax; end
   
   fsteplength = bmax/lsearchsteps;    
   fstepvec = zeros(1,lsearchsteps);
   fstepvec(1:3) = linspace(0,fsteplength,3);
   fstepvec(3:end) = linspace(fstepvec(3),bmax,lsearchsteps-2);
   
   lsearch_go = 1;
   fstep = 0; fsall = 0;
   cnt = 1;  % do not calculate the 0-step (it's already known)
   F = 0;
   F(1) = Fnorm(kk-1);  % 0-step, from the previous iteration
   while lsearch_go & (cnt<lsearchsteps)
     cnt = cnt + 1;
     disp(cnt)
     
     fstep = fstepvec(cnt);
     fsall(cnt) = fstep;
     ss1 = sigma(:,kk-1) + (beta*fstep/bmax).*suunta;
     [F(cnt) U] = LaskeFnormi(ss1,Uel-epsilon1,sigma1,sigma(:,kk-1),args);     
     Uall(:,cnt) = U;
     
     [ff fI] = sort(F(1:cnt));
     if (cnt>2) & ((fI(1)~=1) & (fI(1)~=cnt)), lsearch_go = 0; end     
     %if (fstep+fsteplength*1.5)>bmax, fsteplength = bmax-fstep; end     
   end
   
   % fit polynomial to the data, or select minimum endpoint
   if (fI(1)==1) | (fI(1)==cnt) 
     lmin = fstepvec(fI(1)); 
   else, 
     p = fI(1);
     p = polyfit(fstepvec(p-1:p+1),F(p-1:p+1),2);
     lmin = -.5*p(2)/p(1);           
   end   

   % step vector
   bstep(:,kk) = beta*lmin/bmax;
