% This script compiles pseudo-globally optimal simulation data saved in a
% directory tree into a structure that is easier to manipulate.

blankSlate % clear the workspace

%%%% User inputs %%%%

pdir = [pwd,'/Data/mf0.5']; % directory containing T*/I* subdirectories

updateBest = false; % if true, the script will run GetBestOutputSym in each data directory (for some set of I and T)

%%%%%%%%%%%%%%%%%%%%

list1 = dir([pdir,'/T*']);
list1 = list1([list1.isdir]);% keep only folders
l = 0;

for i = 1:length(list1)
    cdir = [pdir,'/',list1(i).name];
    list2 = dir([cdir,'/I*']);
    list2 = list2([list2.isdir]);
    for j = 1:length(list2)
        cdir2 = [cdir,'/',list2(j).name];
        l = l+1;
        sn = [cdir2,'/BestResultPseudoGlobal'];
        if updateBest
            disp(['----',list1(i).name,', ',list2(j).name,'----'])
            out = GetBestOutputSym(cdir2,'matsavename','BestResultPseudoGlobal',...
                'usepreviousbest',true,'oknoprevbest',true,...
                'BipedDetect','pitch','updatedresult','finemesh',...
                'prevbestname','BestResultPseudoGlobal.mat','Optimizer','objective');
            
            if out.newsol
                A = load(sn);
                animatesolution4SYM(A.outputBest,sn)
                All_results_plotSYM(A.outputBest);
                if datetime(version('-date')) > datetime('February 24, 2020')
                    exportgraphics(gcf,[sn,'.pdf'],'ContentType','vector')
                else
                    saveas(gcf,[sn,'.png'])
                end
            else
                disp('Solution is unchanged')
            end
        end
        A = load(sn); % reload optimal solution
        auxdata = A.outputBest.result.setup.auxdata;
        Optimal(l).auxdata = auxdata;
        Optimal(l).sol = A; % GPOPS-II output
        Optimal(l).obj = A.outputBest.result.objective;
        Optimal(l).I = auxdata.I*4; % convert to murphy number
        Optimal(l).U = auxdata.Uf; % speed
        Optimal(l).D = auxdata.D; % stride length, normalized to body length
        Optimal(l).T = auxdata.D/auxdata.lmax(1)./auxdata.Uf; % Stride period, normalized to front limb length
    end
end

str = datestr(datetime('now'),'yymmddHHMM'); % append current time to savename, to protect overwriting previous results
savename = [pdir,'/BestGaitSolutions', str,'.mat'];

save(savename,'Optimal','pdir')