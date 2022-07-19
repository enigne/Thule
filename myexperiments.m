clear
close all

today = datestr(date(), 'yyyymmdd');
savePath = [today, 'spinup'];
flowmodel = 'MOLHO';

% spin up on a coarse mesh
%steps = [1:2, 4:6];
%resolution = 20e3;
%T = 10000;
%md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'final time', T);

% reinitialize on a finer mesh
%steps = [7];
%resolution = 5e3;
%md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel);


steps = [8];
resolution = 5e3;
md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', 1000);
