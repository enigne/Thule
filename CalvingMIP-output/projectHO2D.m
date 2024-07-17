clear
close all

% Setting {{{
addpath('./');
addpath('../');
addpath('../PostProcessing/');
projectsettings;

glacier = 'Thule';
saveflag = 1;
projPath = ['/totten_1/chenggong/', glacier, '/'];
filename = 'CalvingMIP_Circ_5km_HO_Exp2'; 
%filename = 'cmip-AWI-ex4-G5000-refined-pulse-linear-normal-100.mat'; 
%}}}
% Loading data {{{
org=organizer('repository', [projPath, 'CalvingMIP-output/'], 'prefix', '', 'steps', 1);
disp(['Load the model from ', filename])
if perform(org, [filename, '_save2d']) % {{{
	md = loadmodel(org, filename);
	% project solutions onto the 2D mesh before extrusion  
	md2 = model();
	md2.mesh.numberofelements = md.mesh.numberofelements2d;
	md2.mesh.numberofvertices = md.mesh.numberofvertices2d;
	md2.mesh.elements = md.mesh.elements2d;
	md2.mesh.x = md.mesh.x2d;
	md2.mesh.y = md.mesh.y2d;
	%project solutions to md2
	for i=1:numel(md.results.TransientSolution)
		md2.results.TransientSolution(i).time = md.results.TransientSolution(i).time;
		md2.results.TransientSolution(i).Vx = project2d(md,md.results.TransientSolution(i).Vx,md.mesh.numberoflayers);
		md2.results.TransientSolution(i).Vy = project2d(md,md.results.TransientSolution(i).Vy,md.mesh.numberoflayers);
		md2.results.TransientSolution(i).Thickness = project2d(md,md.results.TransientSolution(i).Thickness,md.mesh.numberoflayers);
		md2.results.TransientSolution(i).MaskIceLevelset = project2d(md,md.results.TransientSolution(i).MaskIceLevelset,md.mesh.numberoflayers);
		md2.results.TransientSolution(i).MaskOceanLevelset = project2d(md,md.results.TransientSolution(i).MaskOceanLevelset,md.mesh.numberoflayers);
	end
	% save md2
	savemodel(org, md2);
end %}}}
%}}}
