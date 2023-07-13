function varargout=runme(varargin)

	%Check inputs {{{
	if nargout>1
		help runme
		error('runme error message: bad usage');
	end
	%recover options
	options=pairoptions(varargin{:});
	% }}}
	%GET cluster name: 'totten'{{{
	clustername = getfieldvalue(options,'cluster name','totten');
	% }}}
	%GET steps: [1]{{{
	steps = getfieldvalue(options,'steps',[1]);
	% }}}
	%GET flow model: 'SSA'{{{
	flowmodel = getfieldvalue(options,'flow model', 'SSA');
	% }}}
	%GET resolution: 5e3{{{
	resolution = getfieldvalue(options,'resolution', 5e3);
	% }}}
	%GET coarse resolution: 10e3{{{
	coarse_resolution = getfieldvalue(options,'coarse resolution', 10e3);
	% }}}
	%GET final time: 1{{{
	finalTime = getfieldvalue(options,'final time', 1);
	% }}}
	%GET relaxation time: 1{{{
	relaxTime = getfieldvalue(options,'relaxation time', 1);
	% }}}
	%GET domain size L: 1e6{{{
	L = getfieldvalue(options,'L', 1e6);
	% }}}
	%GET savepath: '/' {{{
	savePath = getfieldvalue(options,'savePath', '/');
	% }}}
	%GET stabilization for levelset: 5-SUPG {{{
	levelsetStabilization = getfieldvalue(options,'levelset stabilization', 5);
	% }}}
	%GET reinitialization for levelset: 10 {{{
	levelsetReinit = getfieldvalue(options,'levelset reinitialize', 10);
	% }}}
	%GET jobTime for running on supercomputer: 24 hours{{{
	jobTime = getfieldvalue(options,'jobTime', 24);
	% }}}

	%Load some necessary codes {{{
	% main setting
	glacier = 'Thule';

	projPath = ['/totten_1/chenggong/', glacier, '/'];
	% }}}
	%Cluster parameters{{{
	if strcmpi(clustername,'pfe')
		cluster=pfe('numnodes',1,'time',60,'processor','bro','cpuspernode',28,'queue','devel'); %max time is 120 (2hr) and max cpuspernode is 28 for 'bro'
		cluster=pfe('numnodes',1,'time',60,'processor','bro','cpuspernode',28,'queue','normal');
	elseif strcmpi(clustername,'discovery')
		cluster=discovery('numnodes',1,'cpuspernode',64);
		cluster.time = jobTime;
		waitonlock = 0;
	else
		cluster=generic('name',oshostname(),'np', 40);
		waitonlock = Inf;
	end
	clear clustername
	org=organizer('repository',['./Models'],'prefix',['Model_' glacier '_'],'steps',steps); clear steps;
	fprintf(['\n  ========  ' upper(glacier) '  ========\n\n']);
	%}}}
	%Settings and suffix{{{
	suffix = ['_', num2str(resolution/1000, '%.0f'), 'km'];
	coarse_suffix = ['_', num2str(coarse_resolution/1000, '%.0f'), 'km'];
	%}}}za

	%%%%%% Step 1--5
	if perform(org, ['Mesh', suffix])% {{{
		md = triangle(model(), './Exp/Square.exp', resolution);
		md.miscellaneous.name = [glacier, '_', num2str(resolution/1000), 'km'];

		savemodel(org,md);
	end %}}}
	if perform(org, ['Param', suffix])% {{{
		md=loadmodel(org, ['Mesh', suffix]);

		% this does not matter
		md=setflowequation(md,'SSA','all');

		% parameters for CalvingMIP
		md.constants.g = 9.81;				% m s^-2
		md.materials.rho_ice = 917.;		% kg m^-3
		md.materials.rho_water = 1028.;	% kg m^-3
		md.constants.yts = 31556926;     % equivalent to 365.2422*3600*24;
		md.materials.rheology_B = ((2.9377e-9*1e-9)/md.constants.yts)^(-1/3)*ones(md.mesh.numberofvertices,1);  % unit in CalvingMIP: kPa^-3 a^-1, in ISSM: Pa s^1/3
		md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);

		% SMB
		md.smb.mass_balance = 0.3*ones(md.mesh.numberofvertices,1);			% m a^-1 ice eq
		md.basalforcings.groundedice_melting_rate = 0*ones(md.mesh.numberofvertices,1);
		md.basalforcings.floatingice_melting_rate = 0*ones(md.mesh.numberofvertices,1);

		% friction 
		md.friction = frictionweertman();
		% In CalvingMIP setting: u_b=C*tau_b^m, C= 1e-3 m a^-1 kPa^-3
		% In ISSM: tau_b = C^(-1/m) * u_b^(1/m), convert to (m s^-1 Pa^-3)^0.5, 
		md.friction.C = sqrt((1e-12/md.constants.yts)^(-1/3)) * ones(md.mesh.numberofvertices,1); 
		md.friction.m = 3.*ones(md.mesh.numberofelements,1); % m=3

		% geometry
		R=800e3; Bc=900; Bl=-2000; Ba=1100; rc=0;
		% polar coordinates
		r     = sqrt(md.mesh.x.^2 + md.mesh.y.^2);
		theta = atan2(md.mesh.y,md.mesh.x);
		% B calculation
		l=R - cos(2*theta).*R/2 ;
		a=Bc - (Bc-Bl)*(r-rc ).^2./(R-rc ).^2;
		B=Ba*cos(3*pi*r./l)+a ;
		md.geometry.bed = B;
		minimal_thickness = 10;
		md.masstransport.min_thickness=minimal_thickness;
		md.geometry.base = md.geometry.bed;
		md.geometry.surface = md.geometry.base + minimal_thickness;

		% set floating ice    
		pos = (md.geometry.surface<0);
		md.geometry.surface(pos) = (1-md.materials.rho_ice/md.materials.rho_water)*minimal_thickness;
		md.geometry.base(pos) = -md.materials.rho_ice/md.materials.rho_water*minimal_thickness;
		md.geometry.thickness = md.geometry.surface - md.geometry.base;

		% mask
		md = setmask(md,'','');
		md = sethydrostaticmask(md);
		md.mask.ice_levelset = -1*ones(md.mesh.numberofvertices,1);
		md.mask.ice_levelset((r>750e3)) = +1;
		md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

		savemodel(org,md);
	end%}}}
	if perform(org, ['SetBC_', flowmodel, suffix])% {{{

		md=loadmodel(org, ['Param', suffix]);

		% set flow model
		md=setflowequation(md,flowmodel,'all');

		% boundary conditions
		md.stressbalance.spcvx = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.spcvy = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.spcvz = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
		md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

		%Mass transport BC
		md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);

		if strcmp(flowmodel, 'MOLHO')
			disp('  Set boundary conditions for MOLHO');
			md = SetMOLHOBC(md);
		end

		savemodel(org,md);
	end%}}}
	if perform(org, ['Stressbalance_', flowmodel, suffix]) % {{{

		md=loadmodel(org,['SetBC_', flowmodel, suffix]);

		%Set initial speed otherwise solver blows up
		r     = sqrt(md.mesh.x.^2 + md.mesh.y.^2);
		theta = atan2(md.mesh.y,md.mesh.x);
		md.initialization.vx = r.*cos(theta);
		md.initialization.vy = r.*sin(theta);

		% output for MOLHO
		if strcmp(flowmodel, 'MOLHO')
			md.stressbalance.requested_outputs={'default','VxSurface','VySurface','VxShear','VyShear','VxBase','VyBase'};
		end

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.verbose.convergence =1;
		md.cluster = cluster;
		md=solve(md,'sb');

		%Set as initial vx and vy for later
		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end %}}}
	if perform(org, ['Spinup_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Stressbalance_',flowmodel, suffix]);

		% Set parameters
		md.inversion.iscontrol=0;
		md.settings.output_frequency = 100;
		md.timestepping=timesteppingadaptive();
		md.timestepping.time_step_max=100;
		md.timestepping.time_step_min=0.1;
		md.timestepping.start_time=0;
		md.timestepping.final_time=10000;

		% We set the transient parameters
		md.transient.ismovingfront=0;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
		md.stressbalance.requested_outputs={'default'};

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.cluster = cluster;
		md=solve(md,'tr');

		savemodel(org,md);
	end % }}}

	%%%%%% Step 6--10
	% project spinup solution a finer resolution, the reason to not do refinement directly is because the domain is a circle
	if perform(org, ['Reinitialize_', flowmodel, suffix])% {{{
		md=loadmodel(org, ['SetBC_',flowmodel, suffix]);

		if (coarse_resolution>=10e3)
			md_coarse = loadmodel(org,['Spinup_', flowmodel, coarse_suffix]);
		elseif (coarse_resolution<=1e3)
			md_coarse = loadmodel(org,['Pseudo_Relaxation_', flowmodel, coarse_suffix]);
		else
			md_coarse = loadmodel(org,['Exp3_', flowmodel, coarse_suffix]);
		end

		disp(['  Projecting ', num2str(md_coarse.mesh.numberofvertices), ' nodes to a finer mesh with ', num2str(md.mesh.numberofvertices), ' nodes' ])
		minimal_thickness = md.masstransport.min_thickness;

		md.geometry.surface = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Surface, md.mesh.x, md.mesh.y);
		md.geometry.thickness = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Thickness, md.mesh.x, md.mesh.y);
		md.geometry.base = md.geometry.surface - md.geometry.thickness;
		
		% adjust the base
		pos = (md.geometry.base < md.geometry.bed);
		disp(['  Found ', num2str(sum(pos)), ' nodes from the interpolation, where base<bed'])
		md.geometry.base(pos) = md.geometry.bed(pos);
		md.geometry.surface = md.geometry.base + md.geometry.thickness;

		md=sethydrostaticmask(md);

		% set grounded ice base =bed
		pos1 = (md.mask.ocean_levelset>0);
		md.geometry.base(pos1) = md.geometry.bed(pos1);
		md.geometry.surface = md.geometry.base + md.geometry.thickness;

		md=sethydrostaticmask(md);

		md.initialization.vx = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Vx, md.mesh.x, md.mesh.y);
		md.initialization.vy = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Vy, md.mesh.x, md.mesh.y);
		md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

		if strcmp(flowmodel, 'MOLHO')
			md.stressbalance.requested_outputs={'default','VxSurface','VySurface','VxShear','VyShear','VxBase','VyBase'};
		end

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.verbose.convergence =1;
		md.cluster = cluster;

		md=solve(md,'sb');

		%Set as initial vx and vy for later
		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;
		savemodel(org,md);
	end %}}}
	if perform(org, ['Relaxation_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Reinitialize_',flowmodel, suffix]);

		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;

		% Set parameters
		md.inversion.iscontrol=0;
		md.settings.output_frequency = 100;
		md.timestepping=timesteppingadaptive();
		md.timestepping.time_step_max=100;
		md.timestepping.time_step_min=0.01;
		md.timestepping.start_time=0;
		md.timestepping.final_time=relaxTime;

		% We set the transient parameters
		md.transient.ismovingfront=0;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
		md.stressbalance.requested_outputs={'default'};

		md.settings.waitonlock = waitonlock; % do not wait for complete
		if strcmpi(md.cluster.name, 'totten')
			md.miscellaneous.name = ['Thule_transient', suffix];
		else
			md.miscellaneous.name = [savePath];
		end

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		%md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.settings.solver_residue_threshold = 1e-5;
		md.cluster = cluster;
		md=solve(md,'tr', 'runtimename', false);

		savemodel(org,md);
	end % }}}
	if perform(org, ['Exp3_', flowmodel, suffix]) % {{{

		%md=loadmodel(org, ['Relaxation_',flowmodel, suffix]);
		md=loadmodel(org, ['Reinitialize_',flowmodel, suffix]);

		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;

		% Set parameters
		md.inversion.iscontrol=0;
		md.settings.output_frequency = 100;
		md.timestepping=timesteppingadaptive();
		md.timestepping.time_step_max=100;
		md.timestepping.time_step_min=0.01;
		md.timestepping.start_time=0;
		md.timestepping.final_time=relaxTime;

		% We set the transient parameters
		md.transient.ismovingfront=1;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		% set the calving law
		md.calving=calvingcalvingmip();
		md.calving.experiment = 3;  % c=v
		md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
		md.frontalforcings.ablationrate = zeros(md.mesh.numberofvertices,1);
		md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
		md.levelset.stabilization = 5;
		md.levelset.reinit_frequency = 50;

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation', 'MaskOceanLevelset', 'MaskIceLevelset', 'CalvingCalvingrate', 'CalvingMeltingrate', 'GroundedArea','FloatingArea','IcefrontMassFluxLevelset','GroundinglineMassFlux'};
		md.stressbalance.requested_outputs={'default'};

		md.settings.waitonlock = waitonlock; % do not wait for complete
		if strcmpi(md.cluster.name, 'totten')
			md.miscellaneous.name = ['Thule_transient', suffix];
		else
			md.miscellaneous.name = [savePath];
		end

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		md.settings.solver_residue_threshold = 1e-5;
		md.cluster = cluster;
		md=solve(md,'tr', 'runtimename', false);

		savemodel(org,md);
		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end % }}}
	if perform(org, ['Exp4_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Exp3_',flowmodel, suffix]);

		md.initialization.vx = md.results.TransientSolution(end).Vx;
		md.initialization.vy = md.results.TransientSolution(end).Vy;

		% Set parameters
		md.inversion.iscontrol=0;
		md.timestepping=timestepping();
		md.timestepping.start_time=0;
		md.timestepping.final_time=1000;

		% depend on the resolution, 5km->dt=1, 2km->dt=0.4, 1km->dt=0.2 

		%md.timestepping.time_step=cfl_step(md, md.initialization.vx, md.initialization.vy);
		md.timestepping.time_step = 0.5; %1*resolution/5000;
		md.settings.output_frequency = 2; %5000/resolution;

		% We set the transient parameters
		md.transient.ismovingfront=1;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		% set the calving law
		md.calving=calvingcalvingmip();
		md.calving.experiment = 4;  % c=v
		md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);

		% set Wv = -750*sin(2*pi*t/1000);
		Wt = [0:md.timestepping.time_step:500];
		Wv = -750*sin(2*pi*Wt/1000);
		md.frontalforcings.ablationrate = [ones(md.mesh.numberofvertices,1) * Wv; Wt];
		
		md.levelset.spclevelset = NaN(md.mesh.numberofvertices+1,4);
		pos = find(md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos, 1) = md.mask.ice_levelset(pos);
		md.levelset.spclevelset(pos, 2) = md.mask.ice_levelset(pos);
		% set the ice free region at the initial state to be no ice forever
		pos1 = find(md.mask.ice_levelset>0);
		md.levelset.spclevelset(pos1, 3) = sign(md.mask.ice_levelset(pos1));
		md.levelset.spclevelset(pos1, 4) = sign(md.mask.ice_levelset(pos1));
		md.levelset.spclevelset(end,1:4) = [0,500,500.1,1000];

		md.levelset.stabilization = 5;
		md.levelset.reinit_frequency = 50;

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation', 'MaskOceanLevelset', 'MaskIceLevelset', 'CalvingCalvingrate', 'CalvingMeltingrate', 'GroundedArea','FloatingArea','IcefrontMassFluxLevelset','GroundinglineMassFlux'};
		md.stressbalance.requested_outputs={'default'};

		md.settings.waitonlock = waitonlock; % do not wait for complete
		if strcmpi(md.cluster.name, 'totten')
			md.miscellaneous.name = ['Thule_transient', suffix];
		else
			md.miscellaneous.name = [savePath];
		end

		%solve
		%md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.settings.solver_residue_threshold = 1e-5;
		md.cluster = cluster;
%		md.cluster.np=1;
		md.debug.valgrind = 1;
		md.verbose = verbose('all');
		md=solve(md,'tr', 'runtimename', false);

		savemodel(org,md);
		if ~strcmp(savePath, './')
			system(['mkdir -p ./Models/', savePath]);
			system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
		end
	end % }}}



	if perform(org, ['Pseudo_Relaxation_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Reinitialize_',flowmodel, suffix]);

		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;

		% Set parameters
		md.inversion.iscontrol=0;
		md.settings.output_frequency = 100;
		md.timestepping=timesteppingadaptive();
		md.timestepping.time_step_max=1;
		md.timestepping.time_step_min=0.01;
		md.timestepping.start_time=0;
		md.timestepping.final_time=relaxTime;

		% We set the transient parameters
		md.transient.ismovingfront=0;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
		md.stressbalance.requested_outputs={'default'};

		% one nonlinear iter per time step
		md.stressbalance.reltol = NaN;
		md.stressbalance.abstol = NaN;
		md.stressbalance.maxiter = 1;

		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		%md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.settings.solver_residue_threshold = 1e-5;
		md.cluster = cluster;
		md=solve(md,'tr','runtimename',false);

		savemodel(org,md);
	end % }}}
	if perform(org, ['Download_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Pseudo_Relaxation_',flowmodel, suffix]);
		%solve
		md.cluster = discovery('numnodes', 1, 'cpuspernode', 1);

      md.miscellaneous.name = savePath;
      disp(['Downloadng ', savePath, ' from Discovery'])

      md=loadresultsfromcluster(md,'runtimename', savePath);

      savemodel(org,md);

      if ~strcmp(savePath, './')
         system(['mkdir -p ./Models/', savePath]);
         system(['mv ', projPath, '/Models/Model_', glacier, '_', org.steps(org.currentstep).string, '.mat ', projPath, '/Models/', savePath, '/Model_', glacier, '_Transient.mat']);
      end
	end % }}}
	if perform(org, ['Check_', flowmodel, suffix]) % {{{

		md=loadmodel(org, ['Relaxation_SSA', suffix]);

		% set flow model
		md=setflowequation(md,flowmodel,'all');
		if strcmp(flowmodel, 'MOLHO')
			disp('  Set boundary conditions for MOLHO');
			md = SetMOLHOBC(md);
			md.stressbalance.requested_outputs={'default','VxSurface','VySurface','VxShear','VyShear','VxBase','VyBase'};
		end
		% Set parameters
		md.settings.output_frequency = 1;
		md.timestepping=timesteppingadaptive();
		md.timestepping.time_step_max=1;
		md.timestepping.time_step_min=0.01;
		md.timestepping.start_time=0;
		md.timestepping.final_time=1;

		% We set the transient parameters
		md.transient.ismovingfront=0;
		md.transient.isthermal=0;
		md.transient.isstressbalance=1;
		md.transient.ismasstransport=1;
		md.transient.isgroundingline=1;
		md.groundingline.migration = 'SubelementMigration';

		md.verbose.solution=1;
		md.verbose.convergence=0;
		md.cluster = cluster;
		md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
		md.stressbalance.requested_outputs={'default'};

		md.settings.waitonlock = waitonlock; % do not wait for complete
		md.miscellaneous.name = [savePath];

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions('pc_type', 'gamg');
		md.cluster = cluster;
		md=solve(md,'tr', 'runtimename', false);

		savemodel(org,md);
	end % }}}

	%%%%%% Step 11--15
	if perform(org, ['LoadInterpolant', suffix]) % {{{

		md=loadmodel(org, ['Param', suffix]);

		%Start from Hilmar's steady state results
		load('./DATA/SteadyStateInterpolantsThuleMin10km.mat');
		md.geometry.surface = Fs(md.mesh.x,md.mesh.y);
		md.geometry.base    = Fb(md.mesh.x,md.mesh.y);
		pos = find(abs(md.geometry.base - md.geometry.bed)<15);
		md.geometry.base(pos) = md.geometry.bed(pos);
		md.geometry.thickness = md.geometry.surface - md.geometry.base;
		md=sethydrostaticmask(md);

		pos = find(md.mask.ocean_levelset>0);
		md.geometry.base(pos) = md.geometry.bed(pos);
		md.geometry.thickness = md.geometry.surface - md.geometry.base;

		md.mask.ice_levelset = -1*ones(md.mesh.numberofvertices,1);
		md.mask.ice_levelset(find(r>750e3)) = +1;
		md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

		savemodel(org,md);
	end %}}}

	%%%%%% Step 16--20

	%%%%%% Step 21--25

	%%%%%% Step 26--30

	varargout{1} = md;
	return;
