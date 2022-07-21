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
	%GET flow model: 'MOLHO'{{{
	flowmodel = getfieldvalue(options,'flow model', 'MOLHO');
	% }}}
	%GET load from interpolant: 0 {{{
	loadFromInterpolant = getfieldvalue(options,'load from interpolant', 0);
	% }}}
	%GET resolution: 5e3{{{
	resolution = getfieldvalue(options,'resolution', 5e3);
	% }}}
	%GET resolution: 5e3{{{
	coarse_resolution = getfieldvalue(options,'coarse resolution', 20e3);
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
	%GET jobTime for running on supercomputer: 2 hours{{{
	jobTime = getfieldvalue(options,'jobTime', 2);
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
		cluster=discovery('numnodes',1,'cpuspernode',16);
		cluster.time = jobTime;
		waitonlock = 0;
	elseif strcmpi(clustername,'greenplanet')
		cluster=greenplanet('numnodes',1,'cpuspernode',16,'port',0);
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
		md = roundmesh(model(), L, resolution);
		md.miscellaneous.name = [glacier, '_', num2str(resolution/1000), 'km'];

		savemodel(org,md);
	end %}}}
	if perform(org, ['Param', suffix])% {{{
		md=loadmodel(org, ['Mesh', suffix]);

		% this does not matter
		md=setflowequation(md,'SSA','all');

		% parameters for CalvingMIP
		md.constants.yts = 365.2422*3600*24;
		md.constants.g = 9.81;				% m s^-2
		md.materials.rho_ice = 917.;		% kg m^-3
		md.materials.rho_water = 1030.;	% kg m^-3
		md.materials.rheology_B =  ((2.9377e-9*1e-9)/md.constants.yts)^(-1/3)*ones(md.mesh.numberofvertices,1);
		md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);

		% SMB
		md.smb.mass_balance = 0.3*ones(md.mesh.numberofvertices,1);			% m a^-1 ice eq
		md.basalforcings.groundedice_melting_rate = 0*ones(md.mesh.numberofvertices,1);
		md.basalforcings.floatingice_melting_rate = 0*ones(md.mesh.numberofvertices,1);

		% friction 
		md.friction = frictionweertman();
		md.friction.m = 3.*ones(md.mesh.numberofelements,1); % m=3
		md.friction.C = sqrt((1e-12/md.constants.yts)^(-1/3)) * ones(md.mesh.numberofvertices,1); % (m s^-1 Pa^-3)^0.5

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
		minimal_thickness = md.masstransport.min_thickness;
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

		if loadFromInterpolant
			disp('  Use steady state interpolant from Hilmar for the initial condition');
			md=loadmodel(org, ['LoadInterpolant', suffix]);
		else
			md=loadmodel(org, ['Param', suffix]);
		end

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
		md.timestepping.time_step_max=0.1*(cfl_step(md, md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy));
		md.timestepping.time_step_min=0.1;
		md.timestepping.start_time=0;
		md.timestepping.final_time=finalTime;

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

	%%%%%% Step 6--15
	% project spinup solution a finer resolution, the reason to not do refinement directly is because the domain is a circle
	if perform(org, ['Reinitialize_', flowmodel, suffix])% {{{
		md=loadmodel(org, ['Stressbalance_',flowmodel, suffix]);

		md_coarse = loadmodel(org,['Spinup_', flowmodel, coarse_suffix]);

		disp(['  Projecting ', num2str(md_coarse.mesh.numberofvertices), ' nodes to a finer mesh with ', num2str(md.mesh.numberofvertices), ' nodes' ])
		minimal_thickness = md.masstransport.min_thickness;
		
		md.geometry.surface = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Surface, md.mesh.x, md.mesh.y);
		md.geometry.base = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Base, md.mesh.x, md.mesh.y);
		md.geometry.thickness = md.geometry.surface - md.geometry.base;
		md=sethydrostaticmask(md);

		pos = find(md.mask.ocean_levelset>0);
		md.geometry.base(pos) = md.geometry.bed(pos);
		md.geometry.thickness = md.geometry.surface - md.geometry.base;
		md.geometry.thickness(md.geometry.thickness<minimal_thickness) = minimal_thickness;

		md.initialization.vx = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Vx, md.mesh.x, md.mesh.y);
		md.initialization.vy = InterpFromMeshToMesh2d(md_coarse.mesh.elements, md_coarse.mesh.x, md_coarse.mesh.y, md_coarse.results.TransientSolution(end).Vy, md.mesh.x, md.mesh.y);
		md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

		md=solve(md,'sb');
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
		md.timestepping.time_step_max=0.9*(cfl_step(md, md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy));
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

		%solve
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.cluster = cluster;
		md=solve(md,'tr');

		savemodel(org,md);
	end % }}}
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
