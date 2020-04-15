function out = GetBestOutputSym(ParentDir,varargin)
% This function looks through all directories named "rand_guess*" one level
% below Parentdir and pulls out all GPOPSII output files. It then locates
% (by default) the output file that satisfies three conditions:
%   1. Snopt Output is < 10
%   2. mesh error is < mesh tolerance supplied in GPOPS setup
%   3. The objective is minimal among all GPOPSII outputs in the list
% The function then saves a copy of the output, along with information
% about the parent directory and location of the file, in ParentDir
%
% A new option is to specify what to optimize on.
% Parameter: 'optimizer'
% Values: 'objective'(default),'work','forcerate','work+forcerate'
%
% Note: ParentDir MUST be a global (not relative) path!

p = inputParser;
addRequired(p,'ParentDir',@isstr)
validStr = {'supplied','none'};
addParameter(p,'MeshErrorTol','supplied',@(x) (isscalar(x) && isnumeric(x)) || any(strcmpi(x,validStr)) )
addParameter(p,'DisplayResult',true,@islogical)
addParameter(p,'SaveResults',true,@islogical)
addParameter(p,'Verbose',true,@islogical)
addParameter(p,'Optimizer','work+forcerate',@isstr)
addParameter(p,'MATSaveName','BestResult',@isstr)
addParameter(p,'UseInt',true,@islogical)
addParameter(p,'ConvergeCrit','GPOPS',@isstr) %'SNOPT','GPOPS','SNOPTGPOPS'
addParameter(p,'UsePreviousBest',false,@islogical)
addParameter(p,'GuessRange',Inf,@isnumeric)
addParameter(p,'oknoprevbest',false,@islogical)
addParameter(p,'PrevBestName','BestResult.mat',@isstr)
addParameter(p,'ensurealllimbs',false,@islogical)
validStr = {'none','pitch','ensurealllimbs'};
addParameter(p,'BipedDetect','pitch',@(x) any(strcmpi(x,validStr)))
addParameter(p,'IsRecovery',false,@islogical)
parse(p,ParentDir,varargin{:})

MeshErrorTol = p.Results.MeshErrorTol;
DisplayResult = p.Results.DisplayResult;
SaveResults = p.Results.SaveResults;
verbose = p.Results.Verbose;
optimizer = p.Results.Optimizer;
MATSaveName = p.Results.MATSaveName;
ensureAllLimbs = p.Results.ensurealllimbs;
global UseInt
UseInt = p.Results.UseInt;
cc = p.Results.ConvergeCrit;
PrevBest = p.Results.UsePreviousBest;
GuessRange = p.Results.GuessRange;
oknoprevbest = p.Results.oknoprevbest;
PrevBestName = p.Results.PrevBestName;
BipedDetect = p.Results.BipedDetect;
IsRecovery = p.Results.IsRecovery;

UseSuppliedMeshErrTol = false;
if strcmpi(MeshErrorTol,'supplied')
    UseSuppliedMeshErrTol = true;
elseif strcmpi(MeshErrorTol,'none')
    MeshErrorTol = Inf;
end

if any(strcmpi(optimizer,{'prec','recovery','rec'}))
    UseInt = false; % Recovery metric does not use integrals
end

if UseInt
    optimizer = [optimizer,'INT'];
end

% ensure that ParentDir does not end with '/'
if strcmp(ParentDir(end),'/')
    ParentDir = ParentDir(1:end-1);
end
% Find all folders that start with rand_guess

listing = dir([ParentDir,'/rand_guess*']);

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

for curDir = Folders
    curPath = [ParentDir,'/',curDir{1}];
    listing = dir([curPath,'/sol.mat']);
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
        if ensureAllLimbs || strcmpi(BipedDetect,'ensureAllLimbs')
            Ftol = 1e-13;
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
            if isfield(f,'valid')
                valid = f.valid;
            else
                if IsRecovery
                    valid = verifySYM_Rec(output,false);
                else
                    valid = verifySYM(output,false); % check path constraint violation
                end
            end
            crit(4) = ~any(~valid);
            if verbose
                if ~crit(4)
                    disp('Path violation')
                end
            end
        end
        if  all(crit)
            objective = getopt(output,optimizer);
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

function objective = getopt(output,optimizer)
global UseInt
if UseInt
    integrals = output.result.solution.phase.integral;
end

switch optimizer
    case {'prec','recovery','rec'}
        objective = -GPOPS2PREC(output);
    case {'objective','objectiveINT'}
        objective = output.result.objective;
    case 'work'
        [objective,~,~] = CostComponents(output);
    case 'forcerate'
        [~,objective,~] = CostComponents(output);
    case 'work+forcerate'
        [work,forcerate,~] = CostComponents(output);
        objective = work+(forcerate);
    case 'workINT'
        objective = integrals(1);
    case 'forcerateINT'
        objective = integrals(2);
    case 'work+forcerateINT'
        objective = integrals(1) + integrals(2);
    otherwise
        error('Invalid optimizer. Must be ''objective'', ''work'', ''forcerate'', or ''work+forcerate''')
end

end