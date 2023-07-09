function results = ModelToNetCDF(md, varargin)
	%Create netcdf files for an experiment

	%Check inputs {{{
	if nargout>1
		help runme
		error('runme error message: bad usage');
	end
	%recover options
	options=pairoptions(varargin{:});
	% }}}
	%GET directoryname : './/', directory where the outputs are placed {{{
	directoryname = getfieldvalue(options,'directoryname','.//');
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
	%GET author : '' {{{
	author = getfieldvalue(options, 'author', '');
	% }}}

	%Field variables {{{
	disp(['loading field variables..']);

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		results.Time1 = [0]; % Time1
		results.timegrid=[0];
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		results.Time1 = [0:1:1000]; % Time1
		results.timegrid = [0:100:1000]; % Time100
	end

	xmin=-800000; xmax=800000; dx=5000;
	ymin=-800000; ymax=800000; dy=5000;

	results.gridx = xmin:dx:xmax;
	results.gridy = ymin:dy:ymax;
	sizex = numel(results.gridx);
	sizey = numel(results.gridy);

	% allocate matrices
	results.vxmean = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.vymean = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.thickness= nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.mask = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.base = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.calvingrate = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.icemask = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.oceanmask = nan(sizex,sizey,length(results.timegrid)); %y,x,time

	% get the model solution into a matrix
	x = md.mesh.x;
	y = md.mesh.y;
	index = md.mesh.elements;

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		time = 0; 
		vx = md.results.TransientSolution(end).Vx;
		vy = md.results.TransientSolution(end).Vy;
		thickness = md.results.TransientSolution(end).Thickness;
		base = md.results.TransientSolution(end).Base;
		calvingrate = md.results.TransientSolution(end).CalvingCalvingrate; 
		icemask = md.results.TransientSolution(end).MaskIceLevelset;
		oceanmask = md.results.TransientSolution(end).MaskOceanLevelset;
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		time = [0, md.results.TransientSolution(:).time]; % time start at 0, add initial conditions at the beginning
		vx = [[md.initialization.vx], md.results.TransientSolution(:).Vx];
		vy = [[md.initialization.vy], md.results.TransientSolution(:).Vy];
		thickness = [[md.geometry.thickness], md.results.TransientSolution(:).Thickness];
		base = [[md.geometry.base], md.results.TransientSolution(:).Base];
		calvingrate = [[md.results.TransientSolution(1).CalvingCalvingrate], md.results.TransientSolution(:).CalvingCalvingrate]; % we don't have the iniital for this, so extrapolate from t=1
		icemask = [[md.mask.ice_levelset], md.results.TransientSolution(:).MaskIceLevelset];
		oceanmask = [[md.mask.ocean_levelset], md.results.TransientSolution(:).MaskOceanLevelset];
	end
	% set vel/thk in the open ocean to NaN
	vx(icemask>0) = NaN;
	vy(icemask>0) = NaN;
	thickness(icemask>0) = NaN;
	base(icemask>0) = NaN;

	% get these data on the time100 grid
	for i = 1:length(results.timegrid)
		if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
			tid = 1;
		elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
			[~,tid] = min(abs( results.Time1((i-1)*100+1)-time));
		end
		results.vxmean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vx(:, tid), results.gridx, results.gridy, NaN));
		results.vymean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vy(:, tid), results.gridx, results.gridy, NaN));
		results.thickness(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, thickness(:, tid), results.gridx, results.gridy, NaN));
		results.base(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, base(:, tid), results.gridx, results.gridy, NaN));
		results.calvingrate(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, calvingrate(:, tid), results.gridx, results.gridy, NaN));
		results.icemask(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, icemask(:, tid), results.gridx, results.gridy, NaN));
		results.oceanmask(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, oceanmask(:, tid), results.gridx, results.gridy, NaN));
	end

	% mask: grounded=1, floating=2, open ocean=3
	mask = results.icemask;
	mask(results.icemask<0) = 1;		% grounded and floating
	mask(results.icemask>0) = 3;		% open ocean
	mask((results.oceanmask<0) & (results.icemask<0)) = 2;
	results.mask = convertLevelsetsToCalvingMIPMask(results.icemask, results.oceanmask);

	% bed is static
	results.bed = transpose(InterpFromMeshToGrid(index, x, y, md.geometry.bed, results.gridx, results.gridy, NaN));
	%}}}
	%Scalar variables %{{{
	disp(['computing scalar variables..']);

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		results.groundedarea = md.results.TransientSolution(end).GroundedArea; % m^2
		results.floatingarea = md.results.TransientSolution(end).FloatingArea; % m^2
		results.mass = md.results.TransientSolution(end).IceVolume*md.materials.rho_ice; % kg
		results.massaf = md.results.TransientSolution(end).IceVolumeAboveFloatation*md.materials.rho_ice; % kg
		results.totalflux_calving = md.results.TransientSolution(end).IcefrontMassFluxLevelset*10^12; % Gt/yr -> kg/yr
		results.totalflux_groundingline = md.results.TransientSolution(end).GroundinglineMassFlux*10^12;  %Gt/yr -> kg/yr
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		% TODO : add initial values of these quantities from EXP3
		results.groundedarea = [md.results.TransientSolution(:).GroundedArea]; % m^2
		results.floatingarea = [md.results.TransientSolution(:).FloatingArea]; % m^2
		results.mass = [md.results.TransientSolution(:).IceVolume]*md.materials.rho_ice; % kg
		results.massaf = [md.results.TransientSolution(:).IceVolumeAboveFloatation]*md.materials.rho_ice; % kg
		results.totalflux_calving = [md.results.TransientSolution(:).IcefrontMassFluxLevelset]*10^12; % Gt/yr -> kg/yr
		results.totalflux_groundingline = [md.results.TransientSolution(:).GroundinglineMassFlux]*10^12;  %Gt/yr -> kg/yr
	else
		error('not implemented yet');	
	end
	%}}}
	%Profile variables {{{
	disp(['loading profile variables..']);
	if strcmp(expstr,'EXP3') | strcmp(expstr,'EXP4')
		nameList = {'A', 'B', 'C', 'D'};

		% loof through Caprona
		P = readtable('./Results/Caprona_Profiles.csv');
		suffixname = 'Caprona_';
		disp(['  Projecting solutions onto ' suffixname, ' profiles'])
		for i = 1:numel(nameList)
			pfx = P.([suffixname, 'Profile_', nameList{i}, '_X']);
			pfy = P.([suffixname, 'Profile_', nameList{i}, '_Y']);

			% use a temporary profile variable to hold all the info
			pf = [];
			% distance from the start
			pf.distance = P.([suffixname, 'Profile_', nameList{i}, '_S']);

			% project ice thickness, icemask, vx and vy
			pf.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, pfx, pfy);
			pf.vx = InterpFromMeshToMesh2d(index, x, y, vx, pfx, pfy);
			pf.vy = InterpFromMeshToMesh2d(index, x, y, vy, pfx, pfy);
			pf.icemask = InterpFromMeshToMesh2d(index, x, y, icemask, pfx, pfy);
			pf.oceanmask = InterpFromMeshToMesh2d(index, x, y, oceanmask, pfx, pfy);
			pf.mask = convertLevelsetsToCalvingMIPMask(pf.icemask, pf.oceanmask);

			results.profiles.([suffixname, nameList{i}]) = pf;
		end

		% loof through Halbrane
		Q = readtable('./Results/Halbrane_Profiles.csv');
		suffixname = 'Halbrane_';
		disp(['  Projecting solutions onto ' suffixname, ' profiles'])
		for i = 1:numel(nameList)
			pfx = Q.([suffixname, 'Profile_', nameList{i}, '_X']);
			pfy = Q.([suffixname, 'Profile_', nameList{i}, '_Y']);

			% use a temporary profile variable to hold all the info
			pf = [];
			% distance from the start
			pf.distance = Q.([suffixname, 'Profile_', nameList{i}, '_S']);

			% project ice thickness, icemask, vx and vy
			pf.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, pfx, pfy);
			pf.vx = InterpFromMeshToMesh2d(index, x, y, vx, pfx, pfy);
			pf.vy = InterpFromMeshToMesh2d(index, x, y, vy, pfx, pfy);
			pf.icemask = InterpFromMeshToMesh2d(index, x, y, icemask, pfx, pfy);
			pf.oceanmask = InterpFromMeshToMesh2d(index, x, y, oceanmask, pfx, pfy);
			pf.mask = convertLevelsetsToCalvingMIPMask(pf.icemask, pf.oceanmask);

			results.profiles.([suffixname, nameList{i}]) = pf;
		end
	else
		error('not implemented yet');	
	end
	%}}}
	%Create netCDF{{{
	ExpName=[directoryname '/CalvingMIP_' expstr '_' modelname '_' flowequation '_' institution '.nc'];
	disp(['Create netCDF files for the ', expstr, ' at ', ExpName])

	% Delete old file
	if exist(ExpName,'file') == 2
		delete(ExpName);
	end

	% Time{{{
	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		timeStr = 'Time';
		timeExtraStr = 'Time';
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		timeStr = 'Time1';
		timeExtraStr = 'Time100';

		% additional variable
		nccreate(ExpName, timeStr,'Dimensions',{timeStr numel(results.Time1)})
		ncwrite(ExpName,timeStr,results.Time1)
		ncwriteatt(ExpName, timeStr,'units','a');
	end
	nccreate(ExpName, timeExtraStr, 'Dimensions', {timeExtraStr, numel(results.timegrid)})
	ncwrite(ExpName, timeExtraStr, results.timegrid)
	ncwriteatt(ExpName, timeExtraStr,'units','a');
	%}}}
	% Spatial dimensions {{{
	nccreate(ExpName, 'X', 'Dimensions', {'X' 321})
	nccreate(ExpName, 'Y', 'Dimensions', {'Y' 321})
	ncwrite(ExpName, 'X', results.gridx)
	ncwrite(ExpName, 'Y', results.gridy)
	ncwriteatt(ExpName,'X','units','m')
	ncwriteatt(ExpName,'Y','units','m') 
	%}}}
	% Fields: velocity, thickness, masks{{{
	nccreate(ExpName,'xvelmean','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'yvelmean','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'lithk','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'mask','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)})
	nccreate(ExpName,'calverate','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'topg','Dimensions',{'X' 321 'Y' 321})

	ncwrite(ExpName, 'xvelmean', results.vxmean)
	ncwrite(ExpName, 'yvelmean', results.vymean)
	ncwrite(ExpName, 'lithk', results.thickness)
	ncwrite(ExpName, 'mask', results.mask)
	ncwrite(ExpName, 'topg', results.bed)
	ncwrite(ExpName, 'calverate', results.calvingrate)

	ncwriteatt(ExpName,'xvelmean','units','m/a')
	ncwriteatt(ExpName,'yvelmean','units','m/a')
	ncwriteatt(ExpName,'lithk','units','m');
	ncwriteatt(ExpName,'mask','flag_values','1, 2, 3');
	ncwriteatt(ExpName,'topg','units','m');
	ncwriteatt(ExpName,'calverate','units','m/a');

	ncwriteatt(ExpName,'xvelmean','Standard_name','land_ice_vertical_mean_x_velocity')
	ncwriteatt(ExpName,'yvelmean','Standard_name','land_ice_vertical_mean_y_velocity')
	ncwriteatt(ExpName,'lithk','Standard_name','land_ice_thickness');
	ncwriteatt(ExpName,'mask','flag_meanings','1=grounded ice, 2=floating ice, 3=open ocean');
	ncwriteatt(ExpName,'topg','Standard_name','bedrock_altitude');
	ncwriteatt(ExpName,'calverate','Standard_name','calving_rate');
	%}}}
	% scalar variables {{{
	nccreate(ExpName,'iareafl',		'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'iareagr',		'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'lim',				'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'limnsw',			'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'tendlicalvf',	'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'tendligroundf','Dimensions',{timeStr numel(results.Time1)})

	ncwrite(ExpName,'iareafl',results.floatingarea)
	ncwrite(ExpName,'iareagr',results.groundedarea)
	ncwrite(ExpName,'lim',results.mass)
	ncwrite(ExpName,'limnsw',results.massaf)
	ncwrite(ExpName,'tendlicalvf',results.totalflux_calving)
	ncwrite(ExpName,'tendligroundf', results.totalflux_groundingline)

	ncwriteatt(ExpName,'iareafl','units','m^2')
	ncwriteatt(ExpName,'iareagr','units','m^2')
	ncwriteatt(ExpName,'lim','units','kg')
	ncwriteatt(ExpName,'limnsw','units','kg')
	ncwriteatt(ExpName,'tendlicalvf','units','kg/a');
	ncwriteatt(ExpName,'tendligroundf','units','kg/a');

	ncwriteatt(ExpName,'iareafl','Standard_name','grounded_ice_sheet_area')
	ncwriteatt(ExpName,'iareagr','Standard_name','floating_ice_shelf_area')
	ncwriteatt(ExpName,'lim','Standard_name','land_ice_mass')
	ncwriteatt(ExpName,'limnsw','Standard_name','land_ice_mass_not_displacing_sea_water')
	ncwriteatt(ExpName,'tendlicalvf','Standard_name','tendency_of_land_ice_mass_due_to_calving');
	ncwriteatt(ExpName,'tendligroundf','Standard_name','tendency_of_grounded_ice_mass');
	%}}}
	return

	% Profiles
	nccreate(ExpName,'lithkCapA','Dimensions',{'Caprona A' numel(lithkCapA)})
	nccreate(ExpName,'sCapA','Dimensions',{'Caprona A' numel(lithkCapA)})
	nccreate(ExpName,'xvelmeanCapA','Dimensions',{'Caprona A' numel(lithkCapA)},'FillValue',nan)
	nccreate(ExpName,'yvelmeanCapA','Dimensions',{'Caprona A' numel(lithkCapA)},'FillValue',nan)
	nccreate(ExpName,'maskCapA','Dimensions',{'Caprona A' numel(lithkCapA)})

	nccreate(ExpName,'lithkHalA','Dimensions',{'Halbrane A' numel(lithkHalA)})
	nccreate(ExpName,'sHalA','Dimensions',{'Halbrane A' numel(lithkHalA)})
	nccreate(ExpName,'xvelmeanHalA','Dimensions',{'Halbrane A' numel(lithkHalA)},'FillValue',nan)
	nccreate(ExpName,'yvelmeanHalA','Dimensions',{'Halbrane A' numel(lithkHalA)},'FillValue',nan)
	nccreate(ExpName,'maskHalA','Dimensions',{'Halbrane A' numel(lithkHalA)})

	ncwrite(ExpName,'lithkCapA',lithkCapA)
	ncwrite(ExpName,'sCapA',sCapA)
	ncwrite(ExpName,'xvelmeanCapA',xvelmeanCapA)
	ncwrite(ExpName,'xvelmeanCapA',yvelmeanCapA)
	ncwrite(ExpName,'maskCapA',maskCapA)

	ncwrite(ExpName,'lithkHalA',lithkHalA)
	ncwrite(ExpName,'sHalA',sHalA)
	ncwrite(ExpName,'xvelmeanHalA',xvelmeanHalA)
	ncwrite(ExpName,'xvelmeanHalA',yvelmeanHalA)
	ncwrite(ExpName,'maskHalA',maskHalA)

	ncwriteatt(ExpName,'lithkCapA','units','m');
	ncwriteatt(ExpName,'sCapA','units','m');
	ncwriteatt(ExpName,'xvelmeanCapA','units','m/a');
	ncwriteatt(ExpName,'yvelmeanCapA','units','m/a');
	ncwriteatt(ExpName,'maskCapA','flag_values','1, 2, 3');

	ncwriteatt(ExpName,'lithkHalA','units','m');
	ncwriteatt(ExpName,'sHalA','units','m');
	ncwriteatt(ExpName,'xvelmeanHalA','units','m/a');
	ncwriteatt(ExpName,'yvelmeanHalA','units','m/a');
	ncwriteatt(ExpName,'maskHalA','flag_values','1, 2, 3');

	ncwriteatt(ExpName,'lithkCapA','Standard_name','land_ice_thickness_along_profile_A');
	ncwriteatt(ExpName,'sCapA','Standard_name','distance_along_profile_A');
	ncwriteatt(ExpName,'xvelmeanCapA','Standard_name','land_ice_vertical_mean_x_velocity_along_profile_A');
	ncwriteatt(ExpName,'yvelmeanCapA','Standard_name','land_ice_vertical_mean_y_velocity_along_profile_A');
	ncwriteatt(ExpName,'maskCapA','flag_meanings','1=grounded ice, 2=floating ice, 3=open ocean');

	ncwriteatt(ExpName,'lithkHalA','Standard_name','land_ice_thickness_along_profile_A');
	ncwriteatt(ExpName,'sHalA','Standard_name','distance_along_profile_A');
	ncwriteatt(ExpName,'xvelmeanHalA','Standard_name','land_ice_vertical_mean_x_velocity_along_profile_A');
	ncwriteatt(ExpName,'yvelmeanHalA','Standard_name','land_ice_vertical_mean_y_velocity_along_profile_A');
	ncwriteatt(ExpName,'maskHalA','flag_meanings','1=grounded ice, 2=floating ice, 3=open ocean');

	ncdisp(ExpName)

	%}}}
	return
	%Create netcdf {{{
	disp(['creating netcdf...']);

	mode = netcdf.getConstant('NETCDF4');
	mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid=netcdf.create([directoryname '/CalvingMIP_' expstr '_' modelname '_' flowequation '_' institution '.nc'],mode);

	%General attributes
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author', author);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model', modelname);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date', date());

	%Define dimensions
	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		ntime_id  = netcdf.defDim(ncid,'Time',length(results.timegrid));
		time_var_id = netcdf.defVar(ncid,'Time','NC_FLOAT',[ntime_id]);
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		ntime_id  = netcdf.defDim(ncid,'Time100',length(results.timegrid));
		time_var_id = netcdf.defVar(ncid,'Time100','NC_FLOAT',[ntime_id]);

		ntime1_id  = netcdf.defDim(ncid,'Time1',length(results.Time1));
		time1_var_id = netcdf.defVar(ncid,'Time1','NC_FLOAT',[ntime1_id]);
	end
	nx_id  = netcdf.defDim(ncid,'x',length(results.gridx));
	ny_id  = netcdf.defDim(ncid,'y',length(results.gridy));

	disp(['creating netcdf (define variables)...']);
	%Define variables

	xvelmean_var_id= netcdf.defVar(ncid,'xvelmean','NC_FLOAT',[nx_id ny_id ntime_id]);
	netcdf.putAtt(ncid,xvelmean_var_id,'standard_name','land_ice_vertical_mean_x_velocity');
	netcdf.putAtt(ncid,xvelmean_var_id,'units',        'm/a');

	yvelmean_var_id= netcdf.defVar(ncid,'yvelmean','NC_FLOAT',[nx_id ny_id ntime_id]);
	netcdf.putAtt(ncid,yvelmean_var_id,'standard_name','land_ice_vertical_mean_y_velocity');
	netcdf.putAtt(ncid,yvelmean_var_id,'units',        'm/a');

	lithk_var_id= netcdf.defVar(ncid,'lithk','NC_FLOAT',[nx_id ny_id ntime_id]);
	netcdf.putAtt(ncid,lithk_var_id,'standard_name','land_ice_thickness');
	netcdf.putAtt(ncid,lithk_var_id,'units',        'm');

	mask_var_id= netcdf.defVar(ncid,'mask','NC_FLOAT',[nx_id ny_id ntime_id]);
	netcdf.putAtt(ncid,mask_var_id,'standard_name','');
	netcdf.putAtt(ncid,mask_var_id,'units',        'none');

	topg_var_id= netcdf.defVar(ncid,'topg','NC_FLOAT',[nx_id ny_id ntime_id]);
	netcdf.putAtt(ncid,topg_var_id,'standard_name','bedrock_altimetry');
	netcdf.putAtt(ncid,topg_var_id,'units',        'm');

	iareagr_var_id = netcdf.defVar(ncid,'iareagr','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,iareagr_var_id,'standard_name','grounded_ice_sheet_area');
	netcdf.putAtt(ncid,iareagr_var_id,'units',        'm^2');

	iareafl_var_id = netcdf.defVar(ncid,'iareafl','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,iareafl_var_id,'standard_name','floating_ice_sheet_area');
	netcdf.putAtt(ncid,iareafl_var_id,'units',        'm^2');

	lim_var_id = netcdf.defVar(ncid,'lim','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,lim_var_id,'standard_name','land_ice_mass');
	netcdf.putAtt(ncid,lim_var_id,'units',        'kg');

	limnsw_var_id = netcdf.defVar(ncid,'limnsw','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,limnsw_var_id,'standard_name','land_ice_mass_not_displacing_sea_water');
	netcdf.putAtt(ncid,limnsw_var_id,'units',        'kg');

	tendlicalvf_var_id = netcdf.defVar(ncid,'tendlicalvf','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,tendlicalvf_var_id,'standard_name','tendency_of_land_ice_mass_due_to_calving');
	netcdf.putAtt(ncid,tendlicalvf_var_id,'units',        'kg/a');

	tendligroundf_var_id = netcdf.defVar(ncid,'tendligroundf','NC_FLOAT',[ntime_id]);
	netcdf.putAtt(ncid,tendligroundf_var_id,'standard_name','tendency_of_grounded_ice_mass');
	netcdf.putAtt(ncid,tendligroundf_var_id,'units',        'kg/a');

	netcdf.endDef(ncid);

	disp(['creating netcdf (write variables)...']);
	%Write variables
	netcdf.putVar(ncid,xvelmean_var_id,results.vxmean);
	netcdf.putVar(ncid,yvelmean_var_id,results.vymean);
	netcdf.putVar(ncid,lithk_var_id,results.thickness);
	netcdf.putVar(ncid,mask_var_id,results.mask);
	netcdf.putVar(ncid,topg_var_id,results.bed);

	netcdf.putVar(ncid,iareagr_var_id,results.groundedarea);
	netcdf.putVar(ncid,iareafl_var_id,results.floatingarea);
	netcdf.putVar(ncid,lim_var_id,results.mass);
	netcdf.putVar(ncid,limnsw_var_id,results.massaf);
	netcdf.putVar(ncid,tendlicalvf_var_id,results.totalflux_calving);
	netcdf.putVar(ncid,tendligroundf_var_id,results.totalflux_groundingline);

	netcdf.putVar(ncid,time_var_id,results.timescalar);

	%Close netcdf
	netcdf.close(ncid)

	% }}}
