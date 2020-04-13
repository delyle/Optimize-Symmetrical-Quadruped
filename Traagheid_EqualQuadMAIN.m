blankSlate % clear the workspace
run add_necessary_paths % can put the necessary paths in your startup folder instead, or set the relevant paths

%%%% USER INPUTS %%%%
lmax = [1 1]; % Length of forelimbs and hindlimbs (resp.) relative to hip-shoulder length
mf = 0.5; % fore-hind center of mass bias

T = [1 1];%[1 1 1 2 2 2 3 3 3]; % Stride Time 
I = [1 1.5];%[0.75 1 1.5 0.75 1 1.5 0.75 1 1.5]; % Murphy number (normalized pitch moment of inertia). Must be the same length as I

n = 2; % maximum number of guesses to use per T / I combination
pdir = pwd; % Path to parent directory. Guesses will be saved here. Use absolute path
UvsTfun = @(T) (T/2.4).^(-1/0.32); % function to convert stride time to speed. Derived from Alexander & Jayes 1983 J. Zool. Lond. 201:135-152
%%%%%%%%%%%%%%%%%%%%% 


auxdata.mf = mf;
auxdata.lmax = lmax;

pdir1 = strrep([pdir,'/mf',num2str(auxdata.mf,2)],'.','p');
tmpdir = [pdir1,'/tmp2'];

if length(T) ~= length(I)
    error('T and I must be the same length')
end

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
    auxdata.Fr = 1/Tnow^2; % A time constant used in the normalization for the optimizer
    FdotcostScaled = FdotCost/sqrt(0.1)*sqrt(auxdata.Fr);
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
            animatesolution4SYM(out,[savedir,'/sol'],[1 1 1])
        end
    end
    if newsol
        cd(cdir)
        pause(0.5)
        GetBestOutputSym(pwd,'optimizer','objective','matsavename','BestResultTotObj','BipedDetect','pitch'); 
        A = load('BestResultTotObj');
        animatesolution4SYM(A.outputBest,'BestResult')
        All_results_plotSYM(A.outputBest);
        if d > datetime('February 24, 2020')
            exportgraphics(gcf,['BestResult.pdf'],'ContentType','vector')
        else
            saveas(gcf,[savedir,'BestResult.png'])
        end
        
        cd(tmpdir)
        pause(0.5)
        
        save([cdir,'/htoc.mat'],'htoc','auxdata','Fdot*')
    end
end