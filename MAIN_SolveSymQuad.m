% This script iterates through some combination of Murphy number (I/mL^2)
% and T' (T*sqrt(g/(2L))) to find the optimal gait for each combination.
% 
% Default values for the model used in the paper are included
% in the USER INPUTS section. In this example, T and I are specified over a
% sparse grid, and n (the number of guesses per T and I combination) is set
% to a low value. The code can run in a reasonable amount of time. For
% producing the data in the paper, the code was run over a much denser
% grid, and 50-150 guesses were used (depending on I and T).
% 
% --- Running in Parallel ---
% The script can be run with multiple instances of matlab in parallel. It
% will automatically detect whether the current guess has been started and,
% if so, will move on to the next guess.
% 
% Because of this feature, note that if the script is terminated early and 
% rerun, it will start at the next guess after the one in which it was 
% terminated (even if GPOPS-II did not complete its operation for that
% guess).
% 
% If the script is terminated early, a temporary directory with a
% (potentially very large) SNOPT text file will remain in pdir. This can be
% safely deleted manually if the script is not running.

blankSlate % clear the workspace

%%%% USER INPUTS %%%%
lmax = [1 1]; % Length of forelimbs and hindlimbs (resp.) relative to hip-shoulder length
mf = 0.5; % fore-hind center of mass bias

T = [1.5 1.5 1.5 2.5 2.5 2.5 3.5 3.5 3.5]; % Stride Time 
I = [0.75 1 1.5 0.75 1 1.5 0.75 1 1.5]; % Murphy number (normalized pitch moment of inertia). Must be the same length as I

n = 10; % maximum number of guesses to use per T / I combination
pdir = [pwd,'/Data/']; % Path to parent directory. Data will be saved here. Use absolute path
UvsTfun = @(T) (T/2.4).^(-1/0.32); % function to convert stride time to speed. Derived from Alexander & Jayes 1983 J. Zool. Lond. 201:135-152
%%%%%%%%%%%%%%%%%%%%% 



%%%%%%%%%%%%%%%%%%%%%
if length(T) ~= length(I)
    error('T and I must be the same length')
end

auxdata.mf = mf;
auxdata.lmax = lmax;

% make a new directory specifying mf
if strcmp(pdir(end),'/') || strcmp(pdir(end),'\') % remove extraneous '/'
    pdir=  pdir(1:end-1);
end
pdir1 = [pdir,'/mf',num2str(auxdata.mf,2)];

% Make the temporary directory. Appended with the current time to avoid
% conflicts from multiple instances 
tmpdir = [pdir1,'/tmp',datestr(now,'HHMMSS')];
mkdir(tmpdir)
cd(tmpdir)

FdotCost = 3e-5; % Set relative cost of force-rate penalty
auxdata.abounds = pi*[-1 1]; % lower / upper bounds on pitch angle
auxdata.meshtolerance = 1e-3;
auxdata.Fdotmax = 500; % Bounds on absolute rate of change of force

j = (1:n)';
sn = 'rand_guess';
sn = strcat(sn,num2str(j));


d = datetime(version('-date')); % Check release date of current matlab version (for saving figures)

for h = 1:length(T)
    Tnow = T(h);
    auxdata.Uf = UvsTfun(Tnow);
    auxdata.D = Tnow*auxdata.Uf; % Stride length
    auxdata.tau = 1/Tnow^2; % A time constant used in the normalization for the optimizer
    FdotcostScaled = FdotCost/sqrt(0.1)*sqrt(auxdata.tau);
    cdir = [pdir1,'/T',num2str(T(h)),'/I',num2str(I(h))];
    auxdata.I = I(h)/4; % Murphy number converted to the normalization used in the optimization
    mkdir(cdir)
    htoc = [];
    newsol = false;
    for i = 1:size(sn,1)
        auxdata.c = [FdotcostScaled,1,1,1e-3,10,100];
        auxdata.meshMaxIter = 2;
        auxdata.snoptIter = 500;
        auxdata.snopttolerance = 1e-6;
        guess_name = strrep(sn(i,:),' ','');
        savedir = [cdir,'/',guess_name];
        [~,msg] = mkdir(savedir);
        if strcmp(msg,'Directory already exists.')
            disp(msg)
            disp('Skipping to next guess')
        else
            newsol = true;
            tic
            clearvars out*
            close all
            out = SymQuadOptCtrl(auxdata,[]);
            out1 = out;
            All_results_plotSYM(out1);drawnow
            
            
            auxdata.c(5:6) = [100 1000];
            auxdata.meshMaxIter = 3;
            auxdata.snoptIter = 1000;
            auxdata.snopttolerance = 1e-6;
            
            out = SymQuadOptCtrl(auxdata,out1);
            out2 = out;
            All_results_plotSYM(out2);drawnow
            if out2.result.maxerror > auxdata.meshtolerance
                auxdata.meshMaxIter = 8;
                auxdata.c(5:6) = [100 1000];
                out = SymQuadOptCtrl(auxdata,out2);
                out3 = out;
                All_results_plotSYM(out3);drawnow
            end
 
            htoc(i) = toc;
            % Save the figure
            if d > datetime('February 24, 2020')
                exportgraphics(gcf,[savedir,'/sol.pdf'],'ContentType','vector')
            else
                saveas(gcf,[savedir,'/sol.png'])
            end
            save([savedir,'/sol.mat'],'out*')
            
            % refine mesh
            if out.result.nlpinfo < 10 && out.result.maxerror < out.result.setup.mesh.tolerance
                disp(['Refining ',sn(i,:)])
                out = updateFineMesh(out,false); % recompute over a finer mesh
                close all
                All_results_plotSYM(out)
                if d > datetime('February 24, 2020')
                    exportgraphics(gcf,[savedir,'/solFineMesh.pdf'],'ContentType','vector')
                else
                    saveas(gcf,[savedir,'/solFineMesh.png'])
                end
                
                animatesolution4SYM(out,[savedir,'/solFineMesh'],[1 1 1],true,0.02,'PMLIMBS')
                save([savedir,'/solFineMesh.mat'],'out*','auxdata')
            end
        end
    end
    if newsol
        cd(cdir)
        pause(0.5)
        GetBestOutputSym(pwd,'optimizer','objective','updatedresult','finemesh',...
            'matsavename','BestResultPseudoGlobal','BipedDetect','pitch'); 
        A = load('BestResultPseudoGlobal');
        animatesolution4SYM(A.outputBest,'BestResultPseudoGlobal')
        All_results_plotSYM(A.outputBest);
        if d > datetime('February 24, 2020')
            exportgraphics(gcf,[savedir,'/BestResultPseudoGlobal.pdf'],'ContentType','vector')
        else
            saveas(gcf,[savedir,'/BestResultPseudoGlobal.png'])
        end
        
        save([cdir,'/htoc.mat'],'htoc','auxdata','Fdot*')
        cd(tmpdir)
        pause(0.5)
    end
end

cd(pdir)
rmdir(tmpdir,'s')