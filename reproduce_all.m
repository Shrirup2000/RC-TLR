function summary = reproduce_all(profile,outRoot)
%REPRODUCE_ALL Run all RC-TLR experiments currently in this repository.
%   reproduce_all() runs the full numerical sweeps and writes fresh results.
%   reproduce_all('quick') runs a small real-data smoke test of each driver.
%   reproduce_all('full','C:\path\to\results') sets the output directory.
%   Drivers execute one at a time, using their own source folders. The
%   resulting summary distinguishes failed runs from unmet error targets.

    if nargin<1 || isempty(profile), profile='full'; end
    profile=validatestring(profile,{'full','quick'});
    root=fileparts(mfilename('fullpath'));
    if nargin<2 || isempty(outRoot)
        outRoot=fullfile(root,'results',profile);
    end
    if ~isfolder(outRoot), mkdir(outRoot); end
    [ok,info]=fileattrib(outRoot);
    assert(ok && isscalar(info) && info.directory,'Invalid output directory.');
    outRoot=info.Name;

    oldDir=pwd; oldPath=path;
    oldVisibility=get(groot,'DefaultFigureVisible');
    cleanup=onCleanup(@()restore_environment(oldDir,oldPath,oldVisibility)); %#ok<NASGU>
    set(groot,'DefaultFigureVisible','off');
    names=["Newton","Newton_CP","ThreeBox","Biomolecular"];
    folders=["newton","newton","three_box","biomolecular"];
    dirs=fullfile(outRoot,["newton","newton_cp","three_box","biomolecular"]);
    count=numel(names);
    summary=table(names(:),repmat("PENDING",count,1),nan(count,1), ...
        nan(count,1),dirs(:),strings(count,1), ...
        'VariableNames',{'Experiment','Status','CompletedRows','WallSeconds', ...
        'ResultFolder','Message'});
    writetable(summary,fullfile(outRoot,'reproduction_summary.csv'));
    fprintf('\nRC-TLR reproducibility run: %s\nOutput: %s\n',profile,outRoot);

    for k=1:count
        stageDir=char(dirs(k));
        started=tic;
        try
            if ~isfolder(stageDir), mkdir(stageDir); end
            % Avoid collisions between similarly named helpers (for
            % example, lagrange_cheb) in different experiment directories.
            path(oldPath);
            moduleDir=fullfile(root,'experiments',char(folders(k)));
            assert(isfolder(moduleDir),'Missing experiment folder: %s',moduleDir);
            cd(moduleDir);
            addpath(moduleDir,'-begin');
            fprintf('\n[%d/%d] %s\n',k,count,names(k));
            switch names(k)
                case "Newton"
                    assert_driver(moduleDir,'rc_tlr_newton_experiment');
                    [nRows,status,msg]=run_newton(moduleDir,stageDir,profile);
                case "Newton_CP"
                    assert_driver(moduleDir,'rc_tlr_1d_experiment');
                    [nRows,status,msg]=run_newton_cp(moduleDir,stageDir,profile);
                case "ThreeBox"
                    assert_driver(moduleDir,'rc_tlr_threebox_experiment');
                    [nRows,status,msg]=run_three_box(stageDir,profile);
                case "Biomolecular"
                    assert_driver(moduleDir,'run_biomol_paper');
                    [nRows,status,msg]=run_biomolecular(moduleDir,stageDir,profile);
            end
            summary.Status(k)=status;
            summary.CompletedRows(k)=nRows;
            summary.Message(k)=msg;
        catch ME
            summary.Status(k)="FAILED";
            summary.Message(k)=string(ME.message);
            write_failure(fullfile(outRoot,[char(names(k)) '_FAILED.txt']),ME);
            warning('RC_TLR:ReproductionStageFailed','%s: %s',names(k),ME.message);
        end
        summary.WallSeconds(k)=toc(started);
        writetable(summary,fullfile(outRoot,'reproduction_summary.csv'));
    end
    save(fullfile(outRoot,'reproduction_summary.mat'),'summary','profile');
    disp(summary);
    if all(summary.Status=="OK")
        fprintf('All requested RC-TLR experiment stages completed.\n');
    else
        warning('RC_TLR:IncompleteReproduction', ...
            'Some experiment stages are incomplete. Inspect reproduction_summary.csv and the failure reports.');
    end
    fprintf('Result directory: %s\n',outRoot);
end

function [nRows,status,msg]=run_newton(moduleDir,outDir,profile)
    dataDir=fullfile(moduleDir,'data','Newton_CP_-1_1');
    if strcmp(profile,'quick')
        nValues=256; pFull=129; pLong=129;
    else
        nValues=[256 512 1024 2048 4096];
        pFull=[129 257 513 1025 2049 4097 8193 16385];
        pLong=[129 257 513 1025 2049 4097];
    end
    results=table(); firstPass=table();
    grids=table(nValues(:),repmat("PENDING",numel(nValues),1), ...
        zeros(numel(nValues),1),strings(numel(nValues),1), ...
        'VariableNames',{'n','Status','CompletedRows','Message'});
    writetable(grids,fullfile(outDir,'grid_status.csv'));
    % Each grid is independent; save every completed grid immediately so
    % a failure at a larger size does not discard smaller measurements.
    for i=1:numel(nValues)
        n=nValues(i);
        fprintf('Newton grid n=%d\n',n);
        try
            [R,F]=rc_tlr_newton_experiment(dataDir,n,pFull,pLong);
            assert(height(R)==numel(pFull)+numel(pLong),'Incomplete Newton grid sweep.');
            writetable(R,fullfile(outDir,sprintf('newton_n%d.csv',n)));
            writetable(F,fullfile(outDir,sprintf('newton_first_passing_p_n%d.csv',n)));
            if isempty(results)
                results=R; firstPass=F;
            else
                results=[results;R]; firstPass=[firstPass;F]; %#ok<AGROW>
            end
            grids.Status(i)="OK"; grids.CompletedRows(i)=height(R);
        catch ME
            grids.Status(i)="FAILED"; grids.Message(i)=string(ME.message);
            write_failure(fullfile(outDir,sprintf('newton_n%d_FAILED.txt',n)),ME);
            warning('RC_TLR:NewtonGridFailed','n=%d: %s',n,ME.message);
        end
        if ~isempty(results)
            writetable(results,fullfile(outDir,'newton_rc_tlr.csv'));
            writetable(firstPass,fullfile(outDir,'newton_first_passing_p.csv'));
            save(fullfile(outDir,'newton_results.mat'),'results','firstPass');
        end
        writetable(grids,fullfile(outDir,'grid_status.csv'));
    end
    nRows=height(results);
    if all(grids.Status=="OK")
        status="OK"; msg="";
    elseif any(grids.Status=="OK")
        status="PARTIAL"; msg="One or more grids failed; see grid_status.csv.";
    else
        status="FAILED"; msg="No Newton grids completed; see grid_status.csv.";
    end
end

function [nRows,status,msg]=run_newton_cp(moduleDir,outDir,profile)
    opts=struct('dataDir',fullfile(moduleDir,'data','Newton_CP_-1_1'));
    if strcmp(profile,'quick'), opts.pValues=129; end
    [results,model,diagnostics]=rc_tlr_1d_experiment(opts);
    writetable(results,fullfile(outDir,'newton_cp_rc_tlr.csv'));
    save(fullfile(outDir,'newton_cp_results.mat'), ...
        'results','model','diagnostics','opts','-v7.3');
    nRows=height(results); status="OK"; msg="";
end

function [nRows,status,msg]=run_three_box(outDir,profile)
    before=findall(groot,'Type','figure');
    figuresCleanup=onCleanup(@()close_new_figures(before)); %#ok<NASGU>
    if strcmp(profile,'quick')
        [results,comparison]=rc_tlr_threebox_experiment(2:3,1e-7,8);
    else
        [results,comparison]=rc_tlr_threebox_experiment();
    end
    writetable(results,fullfile(outDir,'threebox_rc_tlr.csv'));
    writetable(comparison,fullfile(outDir,'threebox_weighting_comparison.csv'));
    save(fullfile(outDir,'threebox_results.mat'),'results','comparison');
    labels={'Relative Frobenius error','Scalars stored','Construction time (s)'};
    files={'error.png','storage.png','build.png'};
    generated=false(1,numel(labels));
    after=findall(groot,'Type','figure');
    for i=1:numel(after)
        f=after(i);
        if any(f==before), continue; end
        j=find(strcmp(f.Name,labels),1);
        if ~isempty(j)
            exportgraphics(f,fullfile(outDir,files{j}),'Resolution',150);
            generated(j)=true;
        end
    end
    nRows=height(results); status="OK"; msg="";
    if ~all(generated)
        status="PARTIAL"; msg="Tables saved; one or more expected figures were absent.";
    end
end

function [nRows,status,msg]=run_biomolecular(moduleDir,outDir,profile)
    opts=struct('outDir',outDir);
    if strcmp(profile,'quick')
        opts.nValues=256; opts.p=17; opts.pValues=17; opts.NValues=500;
        opts.runRankStudy=false; opts.runInterpolation=false;
        opts.includeAlgorithm32=false; opts.runLattice=false;
        opts.plots=false;
    end
    out=run_biomol_paper(fullfile(moduleDir,'data'),opts);
    nRows=nnz(out.Cases.Status=="OK");
    if out.AllRequestedCasesOK && strlength(out.PlotError)==0
        status="OK"; msg="";
    elseif nRows>0
        status="PARTIAL"; msg="See all_cases.csv and the biomolecular failure reports.";
    else
        status="FAILED"; msg="No biomolecular cases completed.";
    end
end

function assert_driver(moduleDir,name)
    expected=fullfile(moduleDir,[name '.m']);
    assert(strcmp(which(name),expected),'Unexpected MATLAB path resolution for %s.',name);
end

function close_new_figures(before)
    after=findall(groot,'Type','figure');
    before=before(isgraphics(before));
    for i=1:numel(after)
        if ~any(after(i)==before), close(after(i)); end
    end
end

function write_failure(file,ME)
    fid=fopen(file,'w');
    if fid~=-1
        cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
        fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
    end
end

function restore_environment(oldDir,oldPath,oldVisibility)
    cd(oldDir); path(oldPath);
    set(groot,'DefaultFigureVisible',oldVisibility);
end
