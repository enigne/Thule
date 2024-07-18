clear
close all

% Setting 
addpath('./');
addpath('../');
addpath('../PostProcessing/');
projectsettings;

glacier = 'Thule';
saveflag = 1;
projPath = ['/totten_1/chenggong/', glacier, '/'];
filename = 'CalvingMIP_Circ_5km_HO_Exp2_save2d'; frontdata = 'EXP2_HO';
%filename = 'cmip-AWI-ex2-G5000-refined-pulse-linear-normal-100_save2d'; frontdata = 'EXP2_AWI';
%filename = 'cmip-AWI-ex4-G5000-refined-pulse-linear-normal-100_save2d'; frontdata = 'EXP4_AWI';

% Loading data 
org=organizer('repository', [projPath, 'CalvingMIP-output/'], 'prefix', '', 'steps', 0);
disp(['Load the model from ', filename])
md = loadmodel(org, filename);
% get the solutions along the profiles
suffixname = 'Circle';
P = readtable([projPath, 'PostProcessing/Results/', suffixname, '_Profiles.csv']);
profiles = project2Profiles(md, P, suffixname);

% process
names = fieldnames(profiles);
pfnames = names(2:end);
time = profiles.time;
Nt = numel(time);
Np = numel(pfnames);

distance = zeros(Np, Nt);
thickness = zeros(Np, Nt);
vel = zeros(Np, Nt);

for j = 1:Np
	front = getFrontFromProfiles(profiles.(pfnames{j}));
	distance(j,:) = front.distance;
	thickness(j,:) = front.thickness;
	vel(j,:) = front.vel;
end
% save 
save([projPath, 'CalvingMIP-output/Results/', frontdata, '.mat'], 'time', 'distance', 'thickness', 'vel')
