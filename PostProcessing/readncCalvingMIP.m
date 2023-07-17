function results = readncCalvingMIP(varargin)
	%Check inputs {{{
	if nargout>1
		help runme
		error('runme error message: bad usage');
	end
	%recover options
	options=pairoptions(varargin{:});
	% }}}
	%GET directoryname : './/', directory where the outputs are placed {{{
	directoryname = getfieldvalue(options,'directoryname','2_5kmResults');
	if exist(directoryname)~=7,
		error(['directory ' directoryname ' does not exist']);
	end
	% }}}
	%GET EXP : 1, EXP1-EXP4 are currently supported{{{
	expnum = getfieldvalue(options,'EXP', 1);
	expstr = ['EXP', num2str(expnum)];
	% }}}
	%GET modelname : 'ISSM' {{{
	modelname = getfieldvalue(options,'modelname', 'ISSM');
	% }}}
	%GET flowequation : 'SSA' {{{
	flowequation = getfieldvalue(options, 'flowequation', 'SSA');
	% }}}
	%GET institution : 'Dartmouth' {{{
	institution = getfieldvalue(options, 'institution', 'Dartmouth');
	% }}}

% Setting {{{
ExpName=[directoryname '/CalvingMIP_' expstr '_' modelname '_' flowequation '_' institution '.nc'];
%}}}
% Load nc data file {{{
% Time{{{
if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
	timeStr = 'Time';
	timeExtraStr = 'Time';
elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
	timeStr = 'Time1';
	timeExtraStr = 'Time100';
	results.Time1 = ncread(ExpName,timeStr);
end
results.timegrid = ncread(ExpName, timeExtraStr);
disp([' Found ', num2str(numel(results.timegrid)), ' points for ', timeExtraStr])
%}}}
% Spatial dimensions {{{
results.gridx = ncread(ExpName, 'X');
results.gridy = ncread(ExpName, 'Y');
disp(['  The spatial dimension is ', num2str(numel(results.gridx)), ' x ', num2str(numel(results.gridy))])
%}}}
	% Fields: velocity, thickness, masks{{{
	results.vxmean = ncread(ExpName, 'xvelmean');
	results.vymean = ncread(ExpName, 'yvelmean');
	results.thickness = ncread(ExpName, 'lithk');
	results.mask = ncread(ExpName, 'mask');
	results.bed = ncread(ExpName, 'topg');
	results.calvingrate = ncread(ExpName, 'calverate');
	%}}}
	% Scalar variables {{{
	results.floatingarea = ncread(ExpName,'iareafl');
	results.groundedarea = ncread(ExpName,'iareagr');
	results.mass = ncread(ExpName,'lim');
	results.massaf = ncread(ExpName,'limnsw');
	results.totalflux_calving = ncread(ExpName,'tendlicalvf');
	results.totalflux_groundingline = ncread(ExpName,'tendligroundf');
	%}}}
	% Profiles {{{
	nameList = {'A', 'B', 'C', 'D'};
	if strcmp(expstr,'EXP3') | strcmp(expstr,'EXP4')
		fullPre = {'Caprona', 'Halbrane'};		% prefix of the full profile name
		shortPre = {'Cap', 'Hal'};						% prefix of the short profile name
	else
		error('Not implemented')
	end

	for p = 1:numel(fullPre)
		for i = 1:numel(nameList)
			% short and long names for the netCDF file
			sN = [shortPre{p}, nameList{i}];
			fN = [fullPre{p}, ' ', nameList{i}];
			% read varibles from nc files
			pf.thickness = ncread(ExpName, ['lithk', sN]);
			pf.distance = ncread(ExpName, ['s', sN]);
			pf.vx = ncread(ExpName, ['xvelmean', sN]);
			pf.vy = ncread(ExpName, ['yvelmean', sN]);
			pf.mask = ncread(ExpName, ['mask', sN]);
			% Profiles
			results.profiles.([fullPre{p}, '_', nameList{i}]) = pf;
		end
	end
	%}}}
% }}}
