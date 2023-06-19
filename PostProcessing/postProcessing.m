clear
close all
addpath('../');
projectsettings;

glacier = 'Thule';
downloadFromDiscovery = 1;
saveflag = 1;
stepName = 'Transient';
% Setting {{{ 
addpath('./');
projPath = ['/totten_1/chenggong/', glacier, '/'];
steps = [1];
%}}}
% Loading models {{{
[folderList, dataNameList] = getFolderList();
Ndata = length(folderList);
for i = 1:Ndata
	org{i}=organizer('repository', [projPath, 'Models/', folderList{i}], 'prefix', ['Model_' glacier '_'], 'steps', steps);
end

if downloadFromDiscovery  
	% Discovery has limit on max number of connect, so download using serial for-loop
	for i = 1:Ndata
		disp(['---- Downloading the model from Discovery to ', folderList{i}]);
		if perform(org{i}, ['Transient_Discovery_Download'])
			mdList{i} = loadmodel(org{i}, [stepName]);

			mdList{i}.cluster = discovery('numnodes',1,'cpuspernode',1);
			savePath = mdList{i}.miscellaneous.name;

			% download model
			disp(['Downloadng ', savePath, ' from Discovery'])
			mdList{i} = loadresultsfromcluster(mdList{i},'runtimename', savePath);
		end
	end
	parfor i = 1:Ndata
		disp(['Save model to ', savePath])
		% save model
		savemodel(org{i}, mdList{i});
		if ~strcmp(mdList{i}.miscellaneous.name, './')
			system(['mv ', projPath,'/Models/', mdList{i}.miscellaneous.name, '/Model_',glacier,'_Transient_Discovery_Download.mat ', projPath, '/Models/', mdList{i}.miscellaneous.name, '/Model_', glacier, '_', stepName, '.mat']);
		end
	end 
else
	parfor i = 1:Ndata
		disp(['---- Loading the model from ', folderList{i}]);
		mdList{i} = loadmodel(org{i}, [stepName]);
	end
end
%}}}
