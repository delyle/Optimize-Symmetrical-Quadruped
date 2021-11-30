function out = GetBestOutputSym(ParentDir,varargin)
% This function looks through all directories named "rand_guess*" one level
% below Parentdir and pulls out all GPOPSII output files. It then locates
% (by default) the output file that satisfies four conditions:
%   1. Mesh error is < mesh tolerance supplied in GPOPS setup
%   2. SNOPT output flag is < 10, signifying a successful run
%   3. Complementarity constraints are satisfied within tolerance
%   4. The objective is minimal among all GPOPSII outputs in the list
% The function then saves a copy of the output, along with information
% about the parent directory and location of the file, in ParentDir
%
% Note: ParentDir MUST be a global (not relative) path!
%
% --- Name-value pairs ---
% 'MeshErrorTol' -- use tolerances supplied by GPOPS-II, or ignore 
%  'supplied' (default) | 'none' 
%
% 'DiplayResult' -- display path of pseudo-global optimal solution
%  true (default) | false
%
% 'SaveResults' -- save the results as a matlab binary
% true (default) | false
%
% 'Verbose' -- display function progress
% true (default) | false
%
%
% 'Optimizer' -- which cost function to evaluate
% 'objective' (default) | 'work' | 'forcerate' | 'work+forcerate'
% 
% 'objective' uses the GPOPS-II objective. 'work', 'forcerate' and
% 'work+forcerate' calculate those objectives separate from any other
% component(s) of the objective.
%
% 'DirectWork' if true, work is calculated directly from force and length
% changes. If false, slack variables are used.
%
% 'MATSaveName' -- name of .mat file to save pseudo-global optimum
% 'BestResult' (default) | string
%
%
% 'ConvergeCrit' -- whether GPOPS, SNOPT, or both must report success
% 'GPOPS' (default) | 'SNOPT' | 'SNOPTGPOPS'
%
%
% 'GuessRange' -- the range of solutions from individual guesses to
% consider in finding the pseudoglobal optimum. If Inf, then all the
% guesses in ParentDir are considered
% Inf (default) | vector of positive integers
%
%
% 'BipedDetect' -- method to detect (and reject) bipedal solutions
% 'pitch' (default) | 'EnsureAllLimbs' | 'none'
% 
% 'pitch': solution is rejected if the mean absolute pitch angle of the
% body is within pi/32 of pi/2.
%
% 'EnsureAllLimbs': -- Check whether at least one hind limb and one
% forelimb have produced more than 0.01 mg of force in the half cycle
%
%
% 'UsePreviousBest' -- whether to use a previously reported pseudoglobal 
% optimum in ParentDir as the initial standard for other solutions.
% true | false (default)
%
%
% 'PrevBestName' -- filename of the previously found best solution
% 'BestResult' (default) | string
%
%
% 'oknoprevbest' -- whether to set UsePreviousBest to false if the previous 
% best solution cannot be found
% true | false (default)
%
%
% 'RedoConstraintTest' -- whether to check complementarity constraint
% violation, even if the solution is marked as valid
% true | false (default)

p = inputParser;
addRequired(p,'ParentDir',@isstr)
validStr = {'supplied','none'};
addParameter(p,'MeshErrorTol','supplied',@(x) (isscalar(x) && isnumeric(x)) || any(strcmpi(x,validStr)) )
addParameter(p,'DisplayResult',true,@(x) islogical(x) || isscalar(x))
addParameter(p,'SaveResults',true,@(x) islogical(x) || isscalar(x))
addParameter(p,'Verbose',true,@(x) islogical(x) || isscalar(x))
addParameter(p,'Optimizer','work+forcerate',@isstr)
addParameter(p,'DirectWork',false,@(x) islogical(x) || isscalar(x))
addParameter(p,'MATSaveName','BestResult',@isstr)
addParameter(p,'ConvergeCrit','SNOPTGPOPS',@isstr) %'SNOPT','GPOPS','SNOPTGPOPS'
addParameter(p,'GuessRange',Inf,@isnumeric)
validStr = {'none','pitch','ensurealllimbs'};
addParameter(p,'BipedDetect','pitch',@(x) any(strcmpi(x,validStr)))
addParameter(p,'UsePreviousBest',false,@(x) islogical(x) || isscalar(x))
addParameter(p,'PrevBestName','BestResult.mat',@isstr)
addParameter(p,'oknoprevbest',false,@(x) islogical(x) || isscalar(x))
addParameter(p,'redoconstrainttest',false,@(x) islogical(x) || isscalar(x))
addParameter(p,'CostIsNCW',false,@isscalar)
validStr = {'','LRL','LiuRao-Legendre','dFL','dFlmax','finemesh','refinedmesh'};
addParameter(p,'updatedresult','',@(x) any(strcmpi(x,validStr)))
parse(p,ParentDir,varargin{:})

MeshErrorTol = p.Results.MeshErrorTol;
DisplayResult = p.Results.DisplayResult;
SaveResults = p.Results.SaveResults;
verbose = p.Results.Verbose;
optimizer = p.Results.Optimizer;
directWork = p.Results.DirectWork;
MATSaveName = p.Results.MATSaveName;
cc = p.Results.ConvergeCrit;
PrevBest = p.Results.UsePreviousBest;
GuessRange = p.Results.GuessRange;
oknoprevbest = p.Results.oknoprevbest;
PrevBestName = p.Results.PrevBestName;
BipedDetect = p.Results.BipedDetect;
RedoConstraintTest = p.Results.redoconstrainttest;
updatedResult = p.Results.updatedresult;
CostIsNCW = p.Results.CostIsNCW;

UseSuppliedMeshErrTol = false;
if strcmpi(MeshErrorTol,'supplied')
    UseSuppliedMeshErrTol = true;
elseif strcmpi(MeshErrorTol,'none')
    MeshErrorTol = Inf;
end

% ensure that ParentDir does not end with a '/' or '\'
if strcmp(ParentDir(end),filesep)
    ParentDir = ParentDir(1:end-1);
end

% Sanitize ParentDir. Replace all '/' or '\' with the local filesep
ParentDir = strrep(ParentDir,'/',filesep);
ParentDir = strrep(ParentDir,'\',filesep);
% Find all folders that start with rand_guess

listing = dir([ParentDir,filesep,'rand_guess*']);

nGuessesAvailable = length(listing);

if isinf(GuessRange)
    nGuessesUsed = nGuessesAvailable; %#ok<NASGU>
    
    Folders = {listing([listing.isdir]).name};
    
    % Eliminate any folders more than one level below ParentDir
    Folders = Folders(strcmp({listing.folder},ParentDir));
else
    GuessRange = unique(GuessRange); %eliminate duplicates
    nGuessesUsed = length(GuessRange);
    C1 = repmat({'rand_guess'},nGuessesUsed,1);
    C2 = strsplit(num2str(GuessRange(:)'))';
    Folders = strcat(C1,C2)';
end
newsol = false;
if PrevBest
    try
    A = load([ParentDir,'/',PrevBestName]);
    outputBest = A.outputBest;
    objectiveBest = getopt(outputBest,optimizer);
    try
        outputBestFile = A.outputBestFile;
    catch
    end
    outputBestDir = A.outputBestDir;
    catch
       if ~oknoprevbest
           error('No prev best file')
       else
           warning('No prev best file')
           PrevBest = false;
       end
    end
end
if ~PrevBest
    objectiveBest = 1e9;
end
switch upper(updatedResult)
    case {'LRL','LIURAO-LEGENDRE'}
        matname = 'solLRL.mat';
    case {'DFL','DFLMAX'}
        matname = 'soldFlmax.mat';
    case {'REFINEDMESH','FINEMESH'}
        matname = 'solFineMesh.mat';
    otherwise
        matname = 'sol.mat';
end


for curDir = Folders
    curPath = [ParentDir,'/',curDir{1}];
    listing = dir([curPath,'/',matname]);
    if verbose
        disp(['Currently analyzing ',curDir{1}])
    end
    ll = length({listing.name});
    if ll > 1
        warning('More than one outputfile in this directory')
    elseif ll == 0
        warning('Directory is empty')
    else
        i = 1;
        curFile = [listing(i).folder,'/',listing(i).name];
        f = load(curFile);
        output = f.out;
        if UseSuppliedMeshErrTol
            MeshErrorTol = output.result.setup.auxdata.meshtolerance;
        end
        
        crit = logical([1 1 1 1]);
        if contains(upper(cc),'SNOPT')
            crit(1) = output.result.nlpinfo < 10;
            if verbose && ~crit(1)
                disp('SNOPT violation')
            end
        end
        if contains(upper(cc),'GPOPS')
            crit(2) = output.result.maxerror < MeshErrorTol;
            if verbose && ~crit(2)
                disp('MeshTol violation')
            end
        end
        if strcmpi(BipedDetect,'ensureAllLimbs')
            Ftol = 1e-2;
            F = output.result.solution.phase.state(:,7:11);
            a = any(F > Ftol);
            if ~any(a([2,3,5])) || ~any(a([1,4]))
               crit(3) = false;
               if verbose
                   disp('Bipedal: not used')
               end
            end
        end
        if strcmpi(BipedDetect,'pitch')
           theta = output.result.solution.phase.state(:,3);
           if abs(abs(mean(theta)) - pi/2) < pi/32
              if verbose
                  disp('Bipedal (pitch detect): not used')
              end
              crit(3) = false;
           end
        end
        if all(crit)
            if isfield(f,'valid') && ~RedoConstraintTest
                valid = f.valid; % check if this solution is already marked as valid
            elseif CostIsNCW
                valid = verifySYM_Rec(output,false);
            else
                valid = verifySYM(output,false); % check path constraint violation
            end
            crit(4) = ~any(~valid);
            if verbose
                if ~crit(4)
                    disp('Path violation')
                end
            end
        end
        if  all(crit)
            objective = getopt(output,optimizer,directWork);
            if objective < objectiveBest
                newsol = true;
                objectiveBest = objective;
                outputBest = output;
                outputBestFile = curFile;
                outputBestDir = curPath;
            end
        end
    end
end

out.newsol = newsol;
out.objectiveBest = objective;
out.PrevBest = PrevBest;
if SaveResults
    if 1 == exist('outputBest','var')
        save([ParentDir,'/',MATSaveName,'.mat'],'objectiveBest','outputBest','outputBest*','nGuesses*','GuessRange')
        out.outputBest = outputBest;
        out.outputBestDir = outputBestDir;
        out.nGuessesUsed = nGuessesUsed;
        out.nGuessesAvailable = nGuessesAvailable;
    else
        disp('No Feasible Solution')
        save([ParentDir,'/NoFeasibleSolution.mat'],'objectiveBest','nGuesses*','GuessRange')
    end
end

if DisplayResult
    disp(['Best Result is:',char(10),outputBestFile])
end

disp('Process Completed')

end

function objective = getopt(output,optimizer,directWork)

switch optimizer
    case {'prec','recovery','rec'}
        objective = -GPOPS2PREC(output);
    case 'objective'
        objective = output.result.objective;
    case 'work'
        [objective,~,~] = CostComponents(output,directWork);
    case 'forcerate'
        [~,objective,~] = CostComponents(output);
    case 'work+forcerate'
        [work,forcerate,~] = CostComponents(output,directWork);
        objective = work+(forcerate);
    otherwise
        error('Invalid optimizer. Must be ''objective'', ''work'', ''forcerate'', or ''work+forcerate''')
end

end