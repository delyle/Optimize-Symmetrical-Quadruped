function valid = verifySYM_Rec(GPOPSoutput,verbose)

valid = true(1,5);

[F,~,absF,absL,Ldot] = forceveclengthSYM(GPOPSoutput,'grid');

U = GPOPSoutput.result.solution.phase.control;

p = U(:,6);
q = U(:,7);


Ftot = sum(F,3);
v = GPOPSoutput.result.solution.phase.state(:,4:5);
W = dot(Ftot(:,1:2),v,2);

path1 = abs(W - p + q);
path1 = path1(1:end-1); % path of controls sometimes violated at endpoints

crit1 = path1 > 1e-4;
if any(crit1)
   valid(1) = false;
   if verbose
   warning('Violation of W - p + q')
   disp([find(crit1) path1(crit1)])
   end
end

path2 = p.*q;
crit2 = abs(path2) > 1e-4;
if any(crit2(:))
   valid(2) = false;
   if verbose
       warning('Violation of p.*q == 0')
       disp([find(crit2) path2(crit2)])
   end
end

limboff = absF < 0.01;

absLon = absL;
absLon(limboff) = NaN;

aux = GPOPSoutput.result.setup.auxdata;

leg_tol = 0.01;
leg_mult = 1+leg_tol;
crit3 = [absLon(:,[1 4]) > aux.lmax(2)*leg_mult, absLon(:,[2 3 5]) > aux.lmax(1)*leg_mult];
if any(any(crit3))
    valid(3) = false;
    if verbose
    warning('Violation of F*(lmax - l) >= 0')
    disp([find(crit3) absL(crit3)])
    end
end

F_LFt = absF(:,2);
intF_LFl = GPOPSoutput.result.solution.phase.state(:,12);
dF_LFl = GPOPSoutput.result.solution.phase.control(:,3);

path4 = abs(F_LFt.*intF_LFl);
crit4 = path4 > 0.01^2;
if any(crit4)
   valid(4) = false;
   if verbose
   warning('Violation of F_LFt.*intF_LFl')
   disp([find(crit4) path4(crit4)])
   end
end

path5 = abs(F_LFt.*dF_LFl);
path5 = path5(1:end-1);

crit5 = path5 > 0.01;
if any(crit5)
    valid(5) = false;
    if verbose
        warning('Violation of F_LFt.*dF_LFl')
        disp([find(crit5) path4(crit5)])
    end
end

