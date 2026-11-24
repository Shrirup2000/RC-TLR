function test_biomol_rc_tlr()
%TEST_BIOMOL_RC_TLR Synthetic checks, no external data or toolbox needed.
    test_biomol_plan();
    test_rc_tlr_functional();
    rng(17); p=[9 11 13]; R=4;
    U=arrayfun(@(n)randn(n,R),p,'UniformOutput',false);
    CP=struct('lambda',[1;-.3;.4;.2],'U',{U},'weighted',false);
    [wx,wy,wz]=cheb_weight_matrices(p); Gw=weighted_core(CP,wx,wy,wz);
    [I,~,info]=rc_tlr_bounded_nodes(Gw,1e-11,2,.002);
    for d=1:3
        [~,~,~,ordinary]=strong_rrqr_nodes(Gw,d,1e-11,2);
        A=info{d}.sourceMatrix; B=ordinary.sourceMatrix;
        assert(norm(A.'*A-B.'*B,'fro')<1e-11*norm(B.'*B,'fro'));
        assert(info{d}.strongSatisfied && info{d}.accuracySatisfied);
    end
    g=arrayfun(@my_chev2,p,'UniformOutput',false);
    C=selected_core(p,I,CP);
    [X,Y,Z,H]=selected_factor_matrices(g,p,I,info);
    F=rc_tlr_cp_values(CP);
    for k=1:p(3)
        A=reshape(reshape(C,size(C,1)*size(C,2),size(C,3))*Z(k,:).',size(C,1),size(C,2));
        delta=X*A*Y.'-F(:,:,k);
        assert(norm(delta,'fro')<1e-10*norm(F(:,:,k),'fro'));
    end
    for d=1:3, assert(isequal(H{d}(I{d},:),eye(numel(I{d})))); end
    assert(norm(lagrange_cheb(g{1},g{1})-eye(p(1)),'fro')==0);
    root=tempname; mkdir(root);
    fprintf('Synthetic test directory (retained on failure): %s\n',root);
    n=8; g=linspace(-1+1/n,1-1/n,n).';
    LRed=struct('LAM1C',[.7;-.2],'CU1',[ones(n,1),1+g+.3*g.^2], ...
        'CU2',[ones(n,1),1+2*g],'CU3',[ones(n,1),1+.5*g.^3]); %#ok<NASGU>
    Hunif=1; %#ok<NASGU>
    save(fullfile(root,sprintf('n%d_Pn500.mat',n)),'LRed','Hunif');
    opts=struct('nValues',n,'pValues',5,'p',5,'qrTol',1e-11, ...
        'runRankStudy',false,'runInterpolation',false,'includeAlgorithm32',false, ...
        'runLattice',false,'plots',false,'outDir',fullfile(root,'results'));
    out=run_biomol_paper(root,opts);
    assert(out.AllRequestedCasesOK, ...
        'Synthetic integration failed; inspect *_FAILED.txt in %s.',opts.outDir);
    C=out.Cases;
    assert(height(C)==1 && all(C.Method=="RC_TLR"));
    assert(all(C.RelativeError<1e-8));
    assert(all(C.ContinuousStorage==C.CoreEntries+C.p.*(C.rX+C.rY+C.rZ)));
    assert(all(C.GridStorage==C.CoreEntries+C.n.*(C.rX+C.rY+C.rZ)));
    assert(all(C.FactoredGridStorage==C.ContinuousStorage+3*C.n.*C.p));
    assert(C.TotalTime==C.BuildSeconds && C.Repeats==1);
    again=run_biomol_paper(root,opts);
    assert(isequaln(out.Cases,again.Cases));
    assert(isfile(fullfile(opts.outDir,'Table3.csv')));
    assert(isfile(fullfile(opts.outDir,'all_cases.csv')));
    fprintf('RC-TLR synthetic tests passed. Test files: %s\n',root);
end
