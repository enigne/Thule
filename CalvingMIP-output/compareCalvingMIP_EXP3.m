clear 
close all

EXP = 3;

nameList = {'Dartmouth (5km)', 'Dartmouth (1.25km)', 'BAdW (SSA)', 'AWI (HO)'};
fileList = {'../Models/20230705_EXP3_res_5000/Model_Thule_Transient.mat', '../Models/20230718_EXP3_res_1250/Model_Thule_Transient.mat', 'CalvingMIP-Exp3-steady-state-ISSM-BAdW.mat', 'cmip-AWI-ex3-G5000-refined-steady-linear-normal.mat'};

Nf = numel(fileList);
% load models
for i = 1: Nf
	temp = load(fileList{i});
	mdList{i} = temp.md;
end
return

% mesh
figure('Position', [0, 800, 1000, 400])
AXIS = [2e5,4e5,-8e5,-6e5]
for i = 1:Nf
	if i<4
		% SSA
		plotmodel(mdList{i}, 'data', 'mesh', 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS)
	else
		% HO
		plotmodel(mdList{i}, 'data', 'mesh', 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS)
	end
	disp([nameList{i}, '  :  ', num2str((mdList{i}.mesh.numberofelements/1000), '%.0f'), 'k / ', num2str(mdList{i}.mesh.numberofvertices/1000, '%.0f'), 'k'])
end

% steady state
for i = 1:Nf
	if i<4
		% SSA
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(end-1).Vel, 'subplot', [1,Nf+1,i], 'title', nameList{i}, 'caxis', [-1,1])
	else
		% HO
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(end-1).Vel, 'layer', 1, 'subplot', [1,Nf+1,i], 'title', [nameList{i}, ': bottom'], 'caxis', [-1,1])%, 'caxis', [-0.1,0.1])
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(end-1).Vel, 'layer', 7, 'subplot', [1,Nf+1,Nf+1], 'title', [nameList{i}, ': top'], 'caxis', [-1,1])%, 'caxis', [-0.1,0.1])
	end
end


% compare final icemask
figure('Position', [0, 800, 1000, 400])
%AXIS = [4.15e5, 6.5e5,-6.5e5,-4.15e5]
%AXIS = [-8e5, 8e5, -8e5, 8e5]
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
for i = 1:Nf
	if i<4
		% SSA
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset<0, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'caxis', [0,1])
	else
		% HO
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset<0, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'caxis', [0,1])
	end
end

AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
for i = 1:Nf
	if i<4
		% SSA
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'caxis', [-1,1], 'gridded', 1)
	else
		% HO
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'caxis', [-1,1], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
	end
end

% velocity
figure('Position', [0, 800, 1000, 400])
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
AXIS = [-8e5, 8e5, -8e5, 8e5]

for i = 1:Nf
   if i<4
      % SSA
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'caxis',[1000,1200], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   else
      % HO
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'caxis',[1000,1200], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   end
end

% Thickness
figure('Position', [0, 800, 1000, 400])
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
%AXIS = [-8e5, 8e5, -8e5, 8e5]

for i = 1:Nf
   if i<4
      % SSA
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Thickness, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'caxis', [240,350], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   else
      % HO
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Thickness, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'caxis', [240,350], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   end
end
