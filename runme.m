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
	%GET resolution: 5e3{{{
	resolution = getfieldvalue(options,'resolution', 5e3);
	% }}}
	%GET domain size L: 1e6{{{
	L = getfieldvalue(options,'L', 1e6);
	% }}}
	%GET Ice temprature, -8 {{{
	iceTemp = getfieldvalue(options,'iceTemp',-8);
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
		cluster=generic('name',oshostname(),'np', 30);
		waitonlock = Inf;
	end
	clear clustername
	org=organizer('repository',['./Models'],'prefix',['Model_' glacier '_'],'steps',steps); clear steps;
	fprintf(['\n  ========  ' upper(glacier) '  ========\n\n']);
	%}}}
	%set parameters{{{
	%}}}

	%%%%%% Step 1--10
	if perform(org, 'Mesh')% {{{
		md = roundmesh(model(), L, resolution);
		md.miscellaneous.name = [glacier, '_', num2str(resolution/1000), 'km'];

		savemodel(org,md);
	end %}}}
	if perform(org, 'Param')% {{{

		md=loadmodel(org,'Mesh');

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
		minimal_thickness = 10;
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
	if perform(org, 'LoadInterpolant') % {{{

		md=loadmodel(org,'Param');

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
	if perform(org, 'SetBC')% {{{

		md=loadmodel(org,'LoadInterpolant');
		%md=loadmodel(org,'Param');

		%velocity
		md=setflowequation(md,'MOLHO','all');
		%		md=setflowequation(md,'SSA','all');
		md.stressbalance.spcvx = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.spcvy = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.spcvz = NaN(md.mesh.numberofvertices,1);
		md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
		md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

		%Mass transport BC
		md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);

		md = SetMOLHOBC(md);

		savemodel(org,md);
	end%}}}
	if perform(org, 'Stressbalance') % {{{

		md=loadmodel(org,'SetBC');

		%Set initial speed otherwise solver blows up
		r     = sqrt(md.mesh.x.^2 + md.mesh.y.^2);
		theta = atan2(md.mesh.y,md.mesh.x);
		md.initialization.vx = r.*cos(theta);
		md.initialization.vy = r.*sin(theta);
		%md.initialization.vx(:) = 1;
		%md.initialization.vy(:) = 1;

		md.stressbalance.requested_outputs={'default','VxSurface','VySurface','VxShear','VyShear','VxBase','VyBase'};
		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.verbose.convergence =1;
		md.cluster = cluster;
		md=solve(md,'sb');

		%Set as initial vx and vy for later
		md.initialization.vx = md.results.StressbalanceSolution.Vx;
		md.initialization.vy = md.results.StressbalanceSolution.Vy;

		savemodel(org,md);
	end %}}}

	%%%%%% Step 11--15

	%%%%%% Step 16--20

	%%%%%% Step 21--25

	%%%%%% Step 26--30

	varargout{1} = md;
	return;
