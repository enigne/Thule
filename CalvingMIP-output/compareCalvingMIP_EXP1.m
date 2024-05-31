clear 
close all

EXP = 1;

nameList = {'Dartmouth(SSA)', 'Dartmouth(HO)', 'BAdW (SSA)', 'AWI (HO)'};
%fileList = {'CalvingMIP_Exp1_SSA_compressed.mat', 'CalvingMIP_Exp1_HO_5km.mat', 'CalvingMIP-Exp1-steady-state-ISSM-BAdW.mat', 'cmip-AWI-ex1-G5000-refined-relax.mat'};
fileList = {'CalvingMIP_Exp1_SSA_compressed.mat', 'CalvingMIP_Exp1_HO_5km.mat', 'CalvingMIP-Exp1-steady-state-ISSM-BAdW.mat', 'cmip-AWI-ex1-G5000-refined-steady-linear-normal.mat'};
%fileList = {'../Models/20230705_EXP3_res_5000/Model_Thule_Transient.mat', 'CalvingMIP-Exp3-steady-state-ISSM-BAdW.mat', 'cmip-AWI-ex3-G5000-refined-relax.mat'};

Nf = numel(fileList);
% load models
for i = 1: Nf
	temp = load(fileList{i});
	mdList{i} = temp.md;
end
mdList{1}.results.TransientSolution = mdList{1}.results.TransientSolution2;

% mesh
figure('Position', [0, 800, 1000, 400])
AXIS = [2e5,4e5,-8e5,-6e5]
for i = 1:Nf
	if mod(i,2) 
		% SSA
		plotmodel(mdList{i}, 'data', 'mesh', 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS)
	else
		% HO
		plotmodel(mdList{i}, 'data', 'mesh', 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS)
	end
	disp([nameList{i}, '  :  ', num2str((mdList{i}.mesh.numberofelements/1000), '%.0f'), 'k / ', num2str(mdList{i}.mesh.numberofvertices/1000, '%.0f'), 'k'])
end

% equations

% friction law

% ice mask

% calving
for i = 1:2
   if mod(i,2)
      % SSA
      plotmodel(mdList{i}, 'data', mdList{i}.calving.calvingrate, 'subplot', [1,Nf,i], 'title', nameList{i})
   else
      % HO
      plotmodel(mdList{i}, 'data', mdList{i}.calving.calvingrate, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i})
   end
end

% steady state
i = 1;
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution2.Vel, 'subplot', [1,Nf,i], 'title', nameList{i})
i = 2;
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(1).Vel, 'layer', 1, 'subplot', [1,2,1], 'title', 'bottom', 'caxis', [-0.1,0.1])
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(1).Vel, 'layer', 10, 'subplot', [1,2,2], 'title', 'top', 'caxis', [-0.1,0.1])
i = 3;
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(1).Vel, 'subplot', [1,2,1], 'title', '')
i = 4;
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(1).Vel, 'layer', 1, 'subplot', [1,2,1], 'title', 'bottom')%, 'caxis', [-0.1,0.1])
plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel-mdList{i}.results.TransientSolution(1).Vel, 'layer', 7, 'subplot', [1,2,2], 'title', 'top')%, 'caxis', [-0.1,0.1])

% compare final icemask
figure('Position', [0, 800, 1000, 400])
%AXIS = [4.15e5, 6.5e5,-6.5e5,-4.15e5]
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
for i = 1:4
	if mod(i,2) 
		% SSA
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset<0, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS)
	else
		% HO
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).MaskIceLevelset<0, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS)
	end
end


% velocity
figure('Position', [0, 800, 1000, 400])
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
AXIS = [-8e5, 8e5, -8e5, 8e5]

for i = 1:Nf
   if mod(i,2)
      % SSA
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'log', 10, 'caxis',[350,450], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   else
      % HO
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Vel, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'log',10, 'caxis',[350,450], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   end
end

% Thickness
figure('Position', [0, 800, 1000, 400])
AXIS = [-.25e5, 0.25e5,-7.7e5,-7.2e5]
%AXIS = [-8e5, 8e5, -8e5, 8e5]

for i = 1:Nf
   if mod(i,2)
      % SSA
		plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Thickness, 'subplot', [1,Nf,i], 'title', nameList{i}, 'axis', AXIS, 'caxis', [240,300], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   else
      % HO
      plotmodel(mdList{i}, 'data', mdList{i}.results.TransientSolution(end).Thickness, 'subplot', [1,Nf,i], 'layer', 1, 'title', nameList{i}, 'axis', AXIS, 'caxis', [240,300], 'levelset', mdList{i}.results.TransientSolution(end).MaskIceLevelset, 'gridded', 1)
   end
end
