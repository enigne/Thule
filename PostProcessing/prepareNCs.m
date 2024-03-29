clear
close all

%Id = 5;
%outputFolder = './Results/2_5kmResults/';
Id = 4;
outputFolder = './Results/5kmResults/';

% Setting {{{
addpath('./');
addpath('../');
projectsettings;

glacier = 'Thule';
saveflag = 1;
projPath = ['/totten_1/chenggong/', glacier, '/'];
steps = [1];
authorList = 'Cheng, Gong(gong.cheng@dartmouth.edu), Helene Seroussi, Mathieu Morlighem';
%}}}
% Loading models {{{
stepName = 'Transient';
[folderList, dataNameList] = getFolderList(Id);
Ndata = length(folderList);
parfor i = 1:Ndata
	org{i}=organizer('repository', [projPath, 'Models/', folderList{i}], 'prefix', ['Model_' glacier '_'], 'steps', steps);
	if perform(org{i}, 'Transient')
		disp(['---- Loading the model from ', folderList{i}]);
		mdList{i} = loadmodel(org{i}, [stepName]);
	end
end
%}}}
% copy EXP3 final results to EXP4 model if needed {{{
if ~isfield(mdList{2}.results, 'InitialSolution')
	mdList{2}.results.InitialSolution.GroundedArea = mdList{1}.results.TransientSolution(end).GroundedArea;
	mdList{2}.results.InitialSolution.FloatingArea = mdList{1}.results.TransientSolution(end).FloatingArea;
	mdList{2}.results.InitialSolution.IceVolume = mdList{1}.results.TransientSolution(end).IceVolume;
	mdList{2}.results.InitialSolution.IceVolumeAboveFloatation = mdList{1}.results.TransientSolution(end).IceVolumeAboveFloatation;
	mdList{2}.results.InitialSolution.IcefrontMassFluxLevelset = mdList{1}.results.TransientSolution(end).IcefrontMassFluxLevelset;
	mdList{2}.results.InitialSolution.GroundinglineMassFlux = mdList{1}.results.TransientSolution(end).GroundinglineMassFlux;

	savemodel(org{2}, mdList{2});
end
%}}}
% Create NetCDF {{{}}}
results3 = ModelToNetCDF(mdList{1}, 'directoryname', outputFolder, 'EXP', 3, 'author', authorList);
results4 = ModelToNetCDF(mdList{2}, 'directoryname', outputFolder, 'EXP', 4, 'author', authorList);
return
% Generate plots{{{
Cmap = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'};
% set the profiles
NC = 200;
Caprona.A.x = linspace(-390, -590, NC)'*1e3;	Caprona.A.y = linspace(0, 450, NC)'*1e3;
Caprona.B.x = linspace(390, 590, NC)'*1e3;	Caprona.B.y = linspace(0, 450, NC)'*1e3;
Caprona.C.x = linspace(-390, -590, NC)'*1e3; Caprona.C.y = linspace(0, -450, NC)'*1e3;
Caprona.D.x = linspace(390, 590, NC)'*1e3;	Caprona.D.y = linspace(0, -450, NC)'*1e3;
Halbrane.A.x = ones(NC,1)*(-150*1e3);			Halbrane.A.y = linspace(0, 740, NC)'*1e3;
Halbrane.B.x = ones(NC,1)*(150*1e3);			Halbrane.B.y = linspace(0, 740, NC)'*1e3;
Halbrane.C.x = ones(NC,1)*(-150*1e3);			Halbrane.C.y = linspace(0, -740, NC)'*1e3;
Halbrane.D.x = ones(NC,1)*(150*1e3);			Halbrane.D.y = linspace(0, -740, NC)'*1e3;

return
% project solutions to the profile
for nmd=1:Ndata
	% get the model data
	md = mdList{nmd};
	time = [md.results.TransientSolution(:).time];

	% get profiles
	fn = fieldnames(Caprona);
	disp('  Project solutions onto the profiles')
	for i=1:numel(fn)
		Caprona.(fn{i}) = project2Profile(md, Caprona.(fn{i}));
		Halbrane.(fn{i}) = project2Profile(md, Halbrane.(fn{i}));
	end

	% plot ice thickness
	disp('  Start plotting ice thickness and velocity')
	figure('Position', [800, 1000, 1000, 800])
	plotmodel(md, 'data', md.results.TransientSolution(end).Thickness, 'colormap#all', 'parula', ...
		'mask#all', md.results.TransientSolution(end).MaskIceLevelset<0, ...
		'caxis', [0, 1400],...
		'title', 'Ice Thickness (m)', ...
		'subplot', [2,2,1], 'figure', i)
	plotmodel(md, 'data', md.results.TransientSolution(end).Vel, 'colormap#all', 'parula', ...
		'mask#all', md.results.TransientSolution(end).MaskIceLevelset<0, ...
		'caxis', [0, 1200],...
		'title', 'Ice Velocity (m/a)', ...
		'subplot', [2,2,2], 'figure', i)

	% Caprona
	disp('  Plot solutions along the profiles Caprona' )
	subplot(2,2,3)
	i = 1;
	plot(Caprona.(fn{i}).xDist/1e3, Caprona.(fn{i}).bed, 'b')
	hold on
	for i=1:numel(fn)
		plot([Caprona.(fn{i}).xDist; flipud(Caprona.(fn{i}).xDist)]/1e3, [Caprona.(fn{i}).thickness(:,end)+Caprona.(fn{i}).base(:,end); flipud(Caprona.(fn{i}).base(:,end))], 'color', Cmap{i})
	end
	legend(cellfun(@(x) ['Caprona ', x], fn, 'UniformOutput', 0), 'location', 'best')
	xlabel('Distance along the profile (km)')
	ylabel('elevation (m)')
	ylim([-2500, 1000])
	xlim([0, 500])
	title('Caprona ice streams steady state')

	% Halbrane
	disp('  Plot solutions along the profiles Halbrane' )
	subplot(2,2,4)
	i = 1;
	plot(Halbrane.(fn{i}).xDist/1e3, Halbrane.(fn{i}).bed, 'b')
	hold on
	for i=1:numel(fn)
		plot([Halbrane.(fn{i}).xDist; flipud(Halbrane.(fn{i}).xDist)]/1e3, [Halbrane.(fn{i}).base(:,end)+Halbrane.(fn{i}).thickness(:,end); flipud(Halbrane.(fn{i}).base(:,end))], 'color', Cmap{i})
	end
	legend(cellfun(@(x) ['Halbrane ', x], fn, 'UniformOutput', 0))
	xlabel('Distance along the profile (km)')
	ylabel('elevation (m)')
	ylim([-1000, 1500])
	xlim([0, 750])
	title('Halbrane ice streams steady state')
end

%}}}
