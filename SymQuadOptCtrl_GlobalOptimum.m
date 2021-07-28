function SymQuadOptCtrl_GlobalOptimum(ParentDir,auxdata,varargin)
p = inputParser;
addRequired(p,'ParentDir',@ischar)
addRequired(p,'auxdata',@isstruct)

% Values to overwrite auxdata if supplied by user
addParameter(p,'Uf',   [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))
addParameter(p,'D',    [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))
addParameter(p,'mf',   [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))
addParameter(p,'I',    [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))
addParameter(p,'lmax', [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))

addParameter(p,'FdotCost', 3e-5, @(x) isscalar(x) && isnumeric(x))
addParameter(p,'RandIter',   50, @(x) isscalar(x) && isnumeric(x))
addParameter(p,'Recovery',false, @(x) isscalar(x))
addParameter(p,'PointMass',false, @(x) isscalar(x))
addParameter(p,'c5final',[100 1000],@(x) isnumeric(x) && length(x) == 2)
addParameter(p,'forcethirditeration',false, @isscalar)
addParameter(p,'meshmethods',{'hp-PattersonRao','hp-PattersonRao','hp-PattersonRao'},@(x) iscell(x) && length(x) == 3)
parse(p,ParentDir,auxdata,varargin{:})


% Update auxdata if user overwrites values
auxdata = updateAux(auxdata,'Uf',p.Results.Uf);
auxdata = updateAux(auxdata,'D',p.Results.D);
auxdata = updateAux(auxdata,'mf',p.Results.mf);
auxdata = updateAux(auxdata,'I',p.Results.I);
auxdata = updateAux(auxdata,'lmax',p.Results.lmax);



% Retrieve other values
FdotCost = p.Results.FdotCost;
n = p.Results.RandIter;
OptimizeOnRecovery = logical(p.Results.Recovery);
PointMass = logical(p.Results.PointMass);
c5final = p.Results.c5final;
forcethirditeration = logical(p.Results.forcethirditeration);
meshmethods = p.Results.meshmethods;


% Adjust default sigma if user has specified hp-LiuRao-Legendre
if any(strcmp(meshmethods,'hp-LiuRao-Legendre')) && ~isfield(auxdata,'sigma')
    auxdata.sigma = 0.25;
end

auxdata.meshtolerance = 1e-3;

j = (1:n)';
sn = 'rand_guess';
sn = strcat(sn,num2str(j));

auxdata.tau = auxdata.lmax(1)*(auxdata.Uf./auxdata.D).^2; % body length froude number
FdotcostScaled = FdotCost/sqrt(0.1)*sqrt(auxdata.tau);

htoc = [];

if PointMass
    if OptimizeOnRecovery
        GPOPSfun = @SymQuadOptCtrl_pointmass_Rec;
    else
        GPOPSfun = @SymQuadOptCtrl_pointmass;
    end
else
    if OptimizeOnRecovery
        GPOPSfun = @SymQuadOptCtrl_Rec;
    else
        GPOPSfun = @SymQuadOptCtrl;
    end
end

for i = 1:size(sn,1)
    auxdata.c = [FdotcostScaled,1,1,1e-3,c5final*0.1];
    auxdata.meshMaxIter = 2;
    auxdata.snoptIter = 500;
    auxdata.snopttolerance = 1e-6;
    auxdata.meshmethod = meshmethods{1};
    
    guess_name = strrep(sn(i,:),' ',''); % remove whitespace
    savedir = [ParentDir,'/',guess_name];
    [~,msg] = mkdir(savedir);
    if strcmp(msg,'Directory already exists.')
        disp(msg)
        disp('Skipping to next guess')
    else
        tic
        clearvars out*
        close all
        out = GPOPSfun(auxdata,[]);
        if PointMass
            out1 = PointMass2DistMass(out);
        else
            out1 = out;
        end
        All_results_plotSYM(out1);drawnow
        
        
        auxdata.c(5:6) = c5final;
        auxdata.meshMaxIter = 3;
        auxdata.snoptIter = 1000;
        auxdata.snopttolerance = 1e-6;
        auxdata.meshmethod = meshmethods{2};
        
        out = GPOPSfun(auxdata,out);
        if PointMass
            out2 = PointMass2DistMass(out);
        else
            out2 = out;
        end
        All_results_plotSYM(out2);drawnow
        if forcethirditeration || out2.result.maxerror > auxdata.meshtolerance
            auxdata.meshMaxIter = 8;
            auxdata.meshmethod = meshmethods{3};
            auxdata.c(5:6) = c5final;
            out = GPOPSfun(auxdata,out);
            if PointMass
                out3 = PointMass2DistMass(out);
            else
                out3 = out;
            end
            All_results_plotSYM(out3);drawnow
        end
        
        if PointMass
            out = PointMass2DistMass(out);
        end
 
        htoc(i) = toc;
        export_fig([savedir,'/sol.pdf'])
        save([savedir,'/sol.mat'],'out*')
        try
        animatesolution4SYM(out,[savedir,'/sol'])
        catch
        end
    end
end

    cd(ParentDir)
    pause(0.5)
    
    GetBestOutputSym(pwd,'matsavename','BestResultWFR','IsRecovery',OptimizeOnRecovery);
    A = load('BestResultWFR');
    animatesolution4SYM(A.outputBest,'BestResultWFR')
    All_results_plotSYM(A.outputBest);
    export_fig('BestResultWFR.pdf')
    
    GetBestOutputSym(pwd,'optimizer','work','matsavename','BestResultWork','IsRecovery',OptimizeOnRecovery);
    A = load('BestResultWork');
    animatesolution4SYM(A.outputBest,'BestResultWork')
    All_results_plotSYM(A.outputBest);
    export_fig('BestResultWork.pdf')
    
    GetBestOutputSym(pwd,'optimizer','objective','matsavename','BestResultTotObj','IsRecovery',OptimizeOnRecovery);
    A = load('BestResultTotObj');
    animatesolution4SYM(A.outputBest,'BestResultTotObj')
    All_results_plotSYM(A.outputBest);
    export_fig('BestResultTotObj.pdf')
    
    save([ParentDir,'/InputData.mat'],'htoc','auxdata','Fdot*','p')

end

function auxdata = updateAux(auxdata,fieldname,field)

if ~isempty(field)
    auxdata.(fieldname) = field;
end
end