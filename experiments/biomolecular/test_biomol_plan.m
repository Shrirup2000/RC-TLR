function test_biomol_plan()
%TEST_BIOMOL_PLAN Confirm the complete, RC-TLR-only experimental plan.
    [J,c]=biomol_experiment_plan();
    assert(numel(J)==55 && numel(unique({J.id}))==55);
    assert(c.p==129 && isequal(c.nValues,[256 512 1024 2048]));
    assert(isequal(c.NValues,100:100:600));
    assert(isequal(c.pValues,[9 17 33 65 129 257 513 1025]));
    assert(c.qrTol==4e-8 && c.target==1e-7);
    methods={J.method}; systems={J.system};
    selected=strcmp(methods,'RC_TLR');
    spline=strcmp({J.interp},'spline');
    grid=selected & spline & [J.N]==500 & [J.p]==129;
    assert(isequal(sort([J(grid).n]),[256 512 1024 2048]));
    ranks=selected & spline & [J.n]==2048;
    assert(nnz(ranks)==48);
    cases=J(ranks);
    pairs=sortrows([[cases.N].',[cases.p].']);
    [atoms,nodes]=ndgrid(100:100:600,[9 17 33 65 129 257 513 1025]);
    assert(isequal(pairs,sortrows([atoms(:),nodes(:)])));
    storage=ranks & [J.p]==129;
    assert(isequal(sort([J(storage).N]),100:100:600));
    interpolation=selected & [J.n]==2048 & [J.N]==500 & [J.p]==129;
    assert(isequal(sort({J(interpolation).interp}),sort({'spline','linear','nearest'})));
    lattice=selected & strcmp(systems,'lattice');
    assert(nnz(lattice)==1 && J(lattice).n==960 && J(lattice).N==2304);
    full=J(strcmp(methods,'RC_TLR32'));
    assert(numel(full)==1 && full.n==256 && full.N==500 && full.p==129);
    assert(all([J.tol]==c.qrTol));
    assert(all(ismember(methods,{'RC_TLR','RC_TLR32'})));
    short=biomol_experiment_plan(struct('nValues',8,'runRankStudy',false, ...
        'runInterpolation',false,'runLattice',false,'includeAlgorithm32',false,'qrTol',1e-10));
    assert(numel(short)==1 && short.n==8 && short.tol==1e-10);
    fprintf('RC-TLR plan verified: 55 unique cases.\n');
end
