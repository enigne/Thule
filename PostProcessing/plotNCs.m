clear
close all

EXP = 3;
resolutions = [5, 2.5];
for i = 1:numel(resolutions)
	folder = ['./Results/', strrep(num2str(resolutions(i), '%g'),'.','_'), 'kmResults/'];
	names{i} = [num2str(resolutions(i), '%g'), 'km'];
	results{i} = readncCalvingMIP('EXP', EXP, 'directoryname', folder);
end

% Plot transient results{{{
figure('Position', [800, 1000, 1000, 800])
fields = {'groundedarea', 'floatingarea', 'mass', 'massaf', 'totalflux_calving', 'totalflux_groundingline'};
titleList = {'grounded ice sheet area',...
	'floating ice shelf area',...
	'land ice mass',...
	'land ice mass not displacing sea water',...
	'tendency of land ice mass due to calving',...
	'tendency of grounded ice mass'};

for r = 1:numel(resolutions)
	for j = 1:6
		subplot(3,2,j)
		if (EXP == 3) 
			bar(r, results{r}.(fields{j}))
			set(gca, 'XTick', [1:numel(names)])
			set(gca,'XTickLabel', names)
		elseif (EXP ==4) 
			plot(results{r}.Time1, results{r}.(fields{j}))
			%groundedarea)
		else
			error('Not implemented')
		end
		hold on
		title(titleList{j})
	end
end
legend(names, 'location', 'best');
%}}}
% Plot profiles {{{
Cmap = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'};
nameList = {'A', 'B', 'C', 'D'};
fullPre = {'Caprona', 'Halbrane'};     % prefix of the full profile name

for p = 1:numel(fullPre)
	figure('Position', [800, 1000, 1000, 800])
	for i = 1:numel(nameList)
		for r = 1:numel(resolutions)
		% Profiles
		pf = results{r}.profiles.([fullPre{p}, '_', nameList{i}]);
		% short and long names for the netCDF file
		fN = [fullPre{p}, ' ', nameList{i}];
		% write varibles to nc files

		if (EXP == 3) % only show a static figure of the final solution
			subplot(3,1,1)
			plot(pf.distance, pf.thickness, 'color', Cmap{i})
			hold on

			subplot(3,1,2)
			plot(pf.distance, pf.vx, 'color', Cmap{i})
			hold on

			subplot(3,1,3)
			plot(pf.distance, pf.vy, 'color', Cmap{i})
			hold on
		elseif (EXP ==4) % make an animation

		else
			error('Not implemented')
		end
	end
	end
end
%}}}
return
