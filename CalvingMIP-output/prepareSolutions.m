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
EXP = 'EXP4';
%filename = 'CalvingMIP_Circ_5km_HO_Exp2_save2d'; frontdata = 'EXP2_HO';
%filename = 'cmip-AWI-ex2-G5000-refined-pulse-linear-normal-100_save2d'; frontdata = 'EXP2_AWI';
%filename = 'cmip-AWI-ex4-G5000-refined-pulse-linear-normal-100_save2d'; frontdata = 'EXP4_AWI';
%filename = '../Models/20230709_EXP4_res_5000/Model_Thule_Transient'; frontdata = 'EXP4_Dartmouth';
filename = '../Models/20230714_EXP4_res_2500/Model_Thule_Transient'; frontdata = 'EXP4_2_5_Dartmouth';

% Loading data 
org=organizer('repository', [projPath, 'CalvingMIP-output/'], 'prefix', '', 'steps', 0);
disp(['Load the model from ', filename])
md = loadmodel(org, filename);
% get the solutions along the profiles

if strcmp(EXP, 'EXP2')
	suffixname = 'Circle';
	P = readtable([projPath, 'PostProcessing/Results/', suffixname, '_Profiles.csv']);
	profiles = project2Profiles(md, P, suffixname);
elseif strcmp(EXP, 'EXP4')
	suffixname = 'Caprona';
	P = readtable([projPath, 'PostProcessing/Results/', suffixname, '_Profiles.csv']);
	profiles = project2Profiles(md, P, suffixname);
	suffixname = 'Halbrane';
	Q = readtable([projPath, 'PostProcessing/Results/', suffixname, '_Profiles.csv']);
	Q_profiles = project2Profiles(md, Q, suffixname);
	% merge Q_profiles to profiles
	profiles = mergeProfiles(profiles, Q_profiles);
end

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

return

load ../CalvingMIP-output/cmip-AWI-ex4-G5000-refined-pulse-linear-normal-100.mat
mdAwi = md;
load ../Models/20230709_EXP4_res_5000/Model_Thule_Transient.mat
mdDart = md;
mdList = [mdAwi, mdDart];

nameList = {'AWI (HO)', 'Dart (SSA)'};
figure('Position', [0, 800, 500, 400])
for i = 1:numel(nameList)
   md = mdList(i);
   subplot(3,1,1)
   plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).GroundedArea])
   hold on
   subplot(3,1,2)
   plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).FloatingArea])
   hold on
   subplot(3,1,3)
   plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).IceVolume])
   hold on
end
subplot(3,1,1)
xlim([0,1000])
ylabel('Grounded Area')
xlabel('Time (a)')
subplot(3,1,2)
xlim([0,1000])
ylabel('Floating Area')
xlabel('Time (a)')
subplot(3,1,3)
xlim([0,1000])
ylabel('Ice Volume')
xlabel('Time (a)')
legend(nameList,'location', 'best')
set(gcf,'Color','w');


id=3+1;plotmodel(md,'data', md.results.TransientSolution(id).CalvingCalvingrate,'layer#all',7,'mask', md.results.TransientSolution(id).MaskIceLevelset<0,'xlim',[-8e5,8e5],'ylim',[-8e5,8e5],'caxis',ca,'title',['AWI(HO), time=',num2str(md.results.TransientSolution(id).time)])
id=300;plotmodel(mdDart,'data',mdDart.results.TransientSolution(id).CalvingCalvingrate,'mask', mdDart.results.TransientSolution(id).MaskIceLevelset<0,'figure',2,'caxis',ca, 'title',['Dart(HO), time=',num2str(mdDart.results.TransientSolution(id).time)])

ca = [-750,-676];
id=3+1;plotmodel(md,'data',md.results.TransientSolution(id).Vel - md.results.TransientSolution(id).CalvingCalvingrate,'layer#all',7,'mask', md.results.TransientSolution(id).MaskIceLevelset<0,'xlim',[-8e5,8e5],'ylim',[-8e5,8e5],'caxis',ca,'title',['AWI(HO), time=',num2str(md.results.TransientSolution(id).time)])
id=300;plotmodel(mdDart,'data',mdDart.results.TransientSolution(id).Vel - mdDart.results.TransientSolution(id).CalvingCalvingrate,'mask', mdDart.results.TransientSolution(id).MaskIceLevelset<0,'figure',2,'caxis',ca, 'title',['Dart(HO), time=',num2str(mdDart.results.TransientSolution(id).time)])
