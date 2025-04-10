% The main script to run Thule experiments
clear
close all

today = datestr(date(), 'yyyymmdd');

experiments = [13];
flowmodel = 'SSA';

if any(experiments == 1) % exp 1: spin up on a coarse mesh dx=10km {{{
	steps = [1:5];
	resolution = 10e3;
	T = 20000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'final time', T);
end %}}}
if any(experiments == 2) % exp 2: project to 5km mesh and reinitialize {{{
	steps = [1:3, 6];
	resolution = 5e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel);
end %}}}
if any(experiments == 3) % exp 3: relaxation on 5km mesh {{{
	steps = [7];
	resolution = 5e3;
	relaxT = 10000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 4) % exp 4: EXP-3 on 5km mesh {{{
	steps = [8];
	resolution = 5e3;
	relaxT = 10000;
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 20, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 5) % exp 5: EXP-4 on 5km mesh {{{
	steps = [9];
	resolution = 5e3;
	savePath = [today, '_EXP4_res_', num2str(resolution, '%d'), '_NoGLFlux'];
	md = runme('steps', steps,  ...
		'cluster name', 'andes',...
		'jobTime', 5, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel);
end %}}}
if any(experiments == 6) % exp 6: project from 5km to 2.5km mesh and reinitialize {{{
	steps = [1:3, 6];
	coarse_resolution = 5e3;
	resolution = 2.5e3;
	md = runme('steps', steps, ...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 7) % exp 7: EXP-3 on 2.5km mesh {{{
	steps = [8];
	resolution = 2.5e3;
	relaxT = 1000;
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 8) % exp 8: EXP-4 on 2.5km mesh {{{
	steps = [9];
	resolution = 2.5e3;
	savePath = [today, '_EXP4_res_', num2str(resolution, '%d'), '_NoGLFlux'];
	md = runme('steps', steps,  ...
		'cluster name', 'andes',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel);
end %}}}
if any(experiments == 9) % exp 9: project from 2.5km to 1.25km mesh and reinitialize {{{
	steps = [1:3, 6];
	coarse_resolution = 2.5e3;
	resolution = 1.25e3;
	md = runme('steps', steps, ...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 10) % exp 10: EXP-3 on 1.25km mesh {{{
	steps = [8];
	resolution = 1.25e3;
	relaxT = 1000;
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 11) % exp 11: EXP-4 on 1.25km mesh {{{
	steps = [9];
	resolution = 1.25e3;
	savePath = [today, '_EXP4_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel);
end %}}}

if any(experiments == 12) % exp 12: EXP-3 MOLHO on 5km mesh {{{
	steps = [10];
	resolution = 5e3;
	relaxT = 30000;
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'frontera',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 13) % exp 13: EXP-4 MOLHO on 5km mesh {{{
	steps = [11];
	resolution = 5e3;
	savePath = [today, '_EXP4_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'frontera',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution);
end %}}}
return

if any(experiments == 13) % exp 13: project from 5km to 2.5km mesh and reinitialize {{{
	steps = [1:3, 6];
	coarse_resolution = 5e3;
	resolution = 2.5e3;
	md = runme('steps', steps, ...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 14) % exp 14: EXP-3 MOLHO on 2.5km mesh {{{
	steps = [8];
	resolution = 2.5e3;
	relaxT = 1000;
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 40, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 6) % exp 6: relaxation on 2km mesh with moving front (EXP-3) directly {{{
	steps = [8];
	resolution = 2e3;
	relaxT = 1000; % not in used for now
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d'), '_relaxT_', num2str(relaxT, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 10, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}
if any(experiments == 8) % exp 8: relaxation on 1km mesh with moving front (EXP-3) directly {{{
	steps = [8];
	resolution = 1e3;
	relaxT = 1000; % not in used for now
	savePath = [today, '_EXP3_res_', num2str(resolution, '%d'), '_relaxT_', num2str(relaxT, '%d')];
	md = runme('steps', steps,  ...
		'cluster name', 'discovery',...
		'jobTime', 50, ...
		'savePath', [savePath],...
		'resolution', resolution, ...
		'flow model', flowmodel, ...
		'relaxation time', relaxT);
end %}}}


if any(experiments == 6) % exp 6: project from 5km to 2km mesh and reinitialize {{{
	steps = [1:4, 6];
	coarse_resolution = 5e3;
	resolution = 2e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 7) % exp 7: relaxation on 2km mesh {{{
	steps = [7];
	resolution = 2e3;
	relaxT = 5000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 8) % exp 8: relaxation on 2km mesh with maxite=1, discovery {{{
	steps = [8];
	resolution = 2e3;
	relaxT = 10000;
	savePath = [today, '_pseudo_relaxation_2km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 9) % exp 9: project from 2km to 1km mesh and reinitialize {{{
	steps = [1:4, 6];
	coarse_resolution = 2e3;
	resolution = 1e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 10) % exp 10: relaxation on 1km mesh with maxite=1, discovery {{{
	steps = [8];
	resolution = 1e3;
	relaxT = 2000;
	savePath = [today, '_pseudo_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 11) % exp 11: download 1km results from discovery {{{
	steps = [9];
	resolution = 1e3;
	savePath = ['20220808_pseudo_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', 'SSA');
end %}}}
if any(experiments == 12) % exp 12: relaxation on 1km mesh, discovery {{{
	steps = [7];
	resolution = 1e3;
	relaxT = 2000;
	savePath = [today, '_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 13) % exp 13: check not converging on 1km mesh, totten {{{
	steps = [10];
	resolution = 1e3;
%	flowmodel = 'MOLHO';
	savePath = [today, '_check_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'resolution', resolution, 'flow model', flowmodel);
end %}}}

% investigate convergence problem
if any(experiments == 14) % exp 14: set resideu_thres=1e-6 {{{
	steps = [7];
	resolution = 1e3;
	relaxT = 5000;
	savePath = [today, '_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);

	steps = [8];
	savePath = [today, '_pseudo_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 15) % exp 15: download 1km results from discovery {{{
	steps = [9];
	resolution = 1e3;
	savePaths = {
	'20220826_relaxation_1km',...
	};
	for i = 1:length(savePaths)
		md = runme('steps', steps, ...
			'savePath', [savePaths{i}],...
			'cluster name', 'discovery',...
			'resolution', resolution, 'flow model', 'SSA');
	end
end %}}}
