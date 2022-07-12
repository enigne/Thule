steps = [5];

clustername = oshostname();
%Cluster parameters{{{
if strcmpi(clustername,'pfe'),
	interactive = 0;
	cluster=pfe('numnodes',1,'cpuspernode',28,'time',24*5*60,'processor','bro','queue','long','interactive',interactive,'grouplist','s1690');%'g26209');
	cluster=pfe('numnodes',1,'cpuspernode',28,'time',2*60,'processor','bro','queue','devel','interactive',interactive,'grouplist','s1690');
elseif strcmpi(clustername,'greenplanet'),
	cluster=greenplanet('numnodes',3,'cpuspernode',22,'port',0);
elseif strcmpi(clustername,'thwaites') | strcmpi(clustername,'murdo'),
	cluster=generic('name',oshostname(),'np',26);
else
	cluster=generic('name',oshostname(),'np',40);
end
clear clustername interactive
%}}}
org=organizer('repository',['./Models'],'prefix',['Thule_'],'steps',steps); clear steps;

if perform(org,'InitialSetup') % {{{

	resolution = 5e3;
	md=roundmesh(model(),1000e3,resolution);
	md.miscellaneous.name = ['Thule_' num2str(resolution/1000) 'km'];

	%Ice parameters
	md.constants.g = 9.81;
	md.constants.yts = 365.2422*24*3600;
	md.smb.mass_balance = 0.3*ones(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate = 0*ones(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate = 0*ones(md.mesh.numberofvertices,1);
	md.materials.rho_ice = 917.;
	md.materials.rho_water = 1030.;
	md.materials.rheology_B =  ((2.9377*10^-9)*1e3^-3/md.constants.yts)^(-1/3)*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
	md.friction.p=3.*ones(md.mesh.numberofelements,1); % m=3
	md.friction.q=0.*ones(md.mesh.numberofelements,1);
	%md.friction.coefficient = sqrt((3.16*10^6)^-3 *1/(md.constants.yts*1e9)) * ones(md.mesh.numberofvertices,1);
	%md.friction.coefficient = sqrt((3.16*10^6)^-3 *1000*md.constants.yts^3) * ones(md.mesh.numberofvertices,1);
	%md.friction.coefficient = sqrt((3.16*10^6)*1000*md.constants.yts)* ones(md.mesh.numberofvertices,1);
	md.friction.coefficient = sqrt(3.16*10^6)* ones(md.mesh.numberofvertices,1);

	%Set mask and geometry
	% paramters
	R=800e3; Bc=900; Bl=-2000; Ba=1100; rc=0;
	% polar coordinates
	r     = sqrt(md.mesh.x.^2 + md.mesh.y.^2);
	theta = atan2(md.mesh.y,md.mesh.x);
	% B calculation
	l=R - cos(2*theta).*R/2 ;
	a=Bc - (Bc-Bl)*(r-rc ).^2./(R-rc ).^2;
	B=Ba*cos(3*pi*r./l)+a ;
	md.geometry.bed = B;

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

	%velocity
	md=setflowequation(md,'SSA','all');
	md.stressbalance.spcvx = NaN(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy = NaN(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz = NaN(md.mesh.numberofvertices,1);
	md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
	md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

	%Mass transport BC
	md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);

	savemodel(org,md);
end %}}}
if perform(org,'Stressbalance') % {{{

	md=loadmodel(org,'InitialSetup');

	%Set initial speed otherwise solver blows up
	r     = sqrt(md.mesh.x.^2 + md.mesh.y.^2);
	theta = atan2(md.mesh.y,md.mesh.x);
	md.initialization.vx = r.*cos(theta);
	md.initialization.vy = r.*sin(theta);
	%md.initialization.vx(:) = 1;
	%md.initialization.vy(:) = 1;

	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md.verbose.convergence =1;
	md.cluster = cluster;
	md=solve(md,'sb');

	%Set as initial vx and vy for later
	md.initialization.vx = md.results.StressbalanceSolution.Vx;
	md.initialization.vy = md.results.StressbalanceSolution.Vy;

	savemodel(org,md);
end %}}}
if perform(org,'Relaxation') % {{{

	md=loadmodel(org,'Stressbalance');

	% Set parameters
	md.inversion.iscontrol=0;
	md.settings.output_frequency = 50;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=5;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;

	if md.mesh.numberofelements<50e3
		md.timestepping.final_time=5000;
	elseif md.mesh.numberofelements<200e3
		md.timestepping.final_time=200;
	else
		error('not implemented');
	end

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

	plot(cell2mat({md.results.TransientSolution(:).time}),(cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation}) - md.results.TransientSolution(1).IceVolumeAboveFloatation)/md.results.TransientSolution(1).IceVolumeAboveFloatation*100)
	title('Percent change compared to initial state')

	savemodel(org,md);
end % }}}
if perform(org,'Calving=2v') % {{{

	md=loadmodel(org,'Relaxation');

	%Restart model
	md=transientrestart(md);

	% We set the calving model
	md.transient.ismovingfront= 1;
	md.levelset.stabilization = 5; %1 no retreat %5 no retreat %2 
	md.levelset.reinit_frequency = 1;
	md.levelset.migration_max = 2000*365; %2 km / day max
	md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
	pos = find(md.mesh.vertexonboundary);
	md.mask.ice_levelset=reinitializelevelset(md,md.mask.ice_levelset);
	md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);

	%calving
	md.calving=calvingtest();
	md.calving.speedfactor = 2;
	md.settings.solver_residue_threshold = 1e-5;

	%Time
	md.timestepping.time_step_max=1;
	md.timestepping.start_time=0;
	md.timestepping.final_time=200;

	%at least one every 2 years
	deltat = cfl_step(md,md.results.StressbalanceSolution.Vx,md.results.StressbalanceSolution.Vy);
	md.settings.output_frequency = ceil(2/deltat);

	%solve
	md.cluster = cluster;
	md=solve(md,'tr');

	savemodel(org,md);

end % }}}
if perform(org,'Calving=2vMovie') % {{{

	md=loadmodel(org,'Calving=2v');

	%save movie
	plotmodel(md,'data','Vel','icefront','-k','groundingline','-r','caxis',[0 1000],'transient_movie_output',[md.miscellaneous.name '.mp4']);%,...
		%'transient_movie_steps',round(linspace(1,numel(md.results.TransientSolution),150)))
end % }}}
