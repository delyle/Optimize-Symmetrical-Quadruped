function updateResultsFineMesh(pdir,PM_NW,forceredo)
if nargin < 3
    forceredo = false; % redo all refinements in this directory, even if refined solutions already exist
end
PointMass = PM_NW(1);
if length(PM_NW) == 1
    PM_NW(2) = false;
end


list = dir([pdir,'/rand_guess*']);

for i = 1:length(list)
    cdir = [pdir,'/',list(i).name];
    fname = [cdir,'/sol.mat'];
    if forceredo
        skip_guess = false;
    else
        skip_guess = isfile([cdir,'/solFineMesh.mat']);
    end
    if isfile(fname) && ~skip_guess
        A = load(fname);
        input = A.out;
        auxdata = input.result.setup.auxdata;

        if input.result.nlpinfo < 10 && input.result.maxerror < input.result.setup.mesh.tolerance
            disp(['Refining ',list(i).name])
            out = updateFineMesh(input,PM_NW,false);
            if PointMass
                outpm = out; %#ok<NASGU>
                out = PointMass2DistMass(out);
            end
            close all
            All_results_plotSYM(out)
            export_fig([cdir,'/solFineMesh.pdf'])
            if PointMass
               BodyType = 'pointmass';
            else
               BodyType = 'PMLIMBS';
            end
            %animatesolution4SYM(out,[cdir,'/solFineMesh'],[1 1 1],true,0.05,BodyType)
            save([cdir,'/solFineMesh.mat'],'out*','auxdata')
        else
            disp(['SNOPT and/or mesh violation for ',list(i).name])
        end
    else
        disp([list(i).name,' skipped'])
    end
end