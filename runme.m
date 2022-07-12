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
		md=setmask(md,'','');

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



		savemodel(org,md);
	end%}}}

	%%%%%% Step 11--15

	%%%%%% Step 16--20

	%%%%%% Step 21--25

	%%%%%% Step 26--30

	varargout{1} = md;
	return;
