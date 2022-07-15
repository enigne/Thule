clear
close all

today = datestr(date(), 'yyyymmdd');
savePath = [today, 'spinup'];
steps = [1:2, 4:6];
resolution = 20e3;
flowmodel = 'MOLHO';
T = 10000;

md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'final time', T);
