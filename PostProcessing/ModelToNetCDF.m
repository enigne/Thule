function results=ModelToNetCDF(md,directoryname,expnumber,modelname,institution);
	%Create netcdf files for an experiment

	%directoryname: directory where the outputs are placed
	%EXPNUMBER:
	%MODELNAME: ...
	%flowequation
	%INSTITUTION NAME

	if exist(directoryname)~=7,
		error(['directory ' directoryname ' does not exist']);
	end
	
	%Field variables {{{
	disp(['loading field variables..']);

	if strcmp(expnumber,'EXP1') | strcmp(expnumber,'EXP3')
		results.timegrid=[0];
	elseif strcmp(expnumber,'EXP2') | strcmp(expnumber,'EXP4')
		results.timegrid = [0:1:1000];
		results.Time100 = [0:100:1000];
	end

	xmin=-800000; xmax=800000; dx=5000;
	ymin=-800000; ymax=800000; dy=5000;

	results.gridx = xmin:dx:xmax;
	results.gridy = ymin:dy:ymax;
	sizex = numel(results.gridx);
	sizey = numel(results.gridy);
	
	results.vxmean = nan(sizex,sizey,length(results.Time100)); %y,x,time
	results.vymean = nan(sizex,sizey,length(results.Time100)); %y,x,time
	results.thickness= nan(sizex,sizey,length(results.Time100)); %y,x,time
	results.mask= nan(sizex,sizey,length(results.Time100)); %y,x,time
	results.base= nan(sizex,sizey,length(results.Time100)); %y,x,time
	results.calvingrate= nan(sizex,sizey,length(results.Time100)); %y,x,time

	% get the model solution into a matrix
	x = md.mesh.x;
	y = md.mesh.y;
	index = md.mesh.elements;
	time = [md.results.TransientSolution(:).time];
	vx = [md.results.TransientSolution(:).Vx];
	vy = [md.results.TransientSolution(:).Vy];
	thickness = [md.results.TransientSolution(:).Thickness];
	base = [md.results.TransientSolution(:).Base];
	calvingrate = [md.results.TransientSolution(:).CalvingCalvingrate];
	icemask = [md.results.TransientSolution(:).MaskIceLevelset];

	% set vel/thk in the open ocean to NaN
	vx(icemask>0) = NaN;
	vy(icemask>0) = NaN;
	thickness(icemask>0) = NaN;
	base(icemask>0) = NaN;

	% get these data on the time100 grid
	for i = 1:length(results.Time100)
		[~,tid] = min(abs( results.timegrid((i-1)*100+1)-time));
		results.vxmean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vx(:, tid), results.gridx, results.gridy, NaN));
		results.vymean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vy(:, tid), results.gridx, results.gridy, NaN));
		results.thickness(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, thickness(:, tid), results.gridx, results.gridy, NaN));
		results.base(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, base(:, tid), results.gridx, results.gridy, NaN));
		results.calvingrate(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, calvingrate(:, tid), results.gridx, results.gridy, NaN));
	end
	return
	
	[vymean]=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Vy,results.gridx,results.gridy,NaN);
	results.vymean=transpose(vymean);
	
	[thickness]=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).Thickness,results.gridx,results.gridy,NaN);
	results.thickness=transpose(thickness);

	mask=-md.results.TransientSolution(end).MaskIceLevelset;
	mask(find(mask>0))=1;
	mask(find(mask<0))=3;
	pos=find(md.results.TransientSolution(end).MaskOceanLevelset<0 & md.results.TransientSolution(end).MaskIceLevelset<0); mask(pos)=2;

	F = scatteredInterpolant(md.mesh.x,md.mesh.y,mask,'nearest','none');
	[xq,yq]=meshgrid(results.gridx,results.gridy);
	icemask = F(xq,yq);
	%[icemask]=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,mask,results.gridx,results.gridy,NaN);
	results.mask=transpose(icemask);

	[bed]=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed,results.gridx,results.gridy,NaN);
	results.bed=transpose(bed);

	%}}}

	%Scalar variables %{{{
	disp(['loading scalar variables..']);
	
	if strcmp(expnumber,'EXP1'),
		results.timescalar=[0];
	else
		results.timescalar=[1:1:101];
	end
	
	results.groundedarea=[];
	results.floatingndarea=[];
	results.mass=[];
	results.massaf=[];
	results.totalflux_calving=[];
	results.totalflux_groundingline=[];

	if strcmp(expnumber,'EXP1'),
		results.groundedarea=md.results.TransientSolution(end).GroundedArea;
		results.floatingarea=md.results.TransientSolution(end).FloatingArea;
		results.mass=md.results.TransientSolution(end).IceVolume*md.materials.rho_ice;
		results.massaf=md.results.TransientSolution(end).IceVolumeAboveFloatation*md.materials.rho_ice;
		results.totalflux_calving=md.results.TransientSolution(end).IcefrontMassFluxLevelset*10^12; % kg/yr
		results.totalflux_groundingline=md.results.TransientSolution(end).GroundinglineMassFlux*10^12;  %kg/yr
	else
		error('not implemented yet');	
	end
	
	%}}}
	
%	%Profile variables {{{
%	disp(['loading profile variables..']);
%	
%	if strcmp(expnumber,'EXP1'),
%		results.timegrid=[0];
%	else
%		results.timegrid=[1:5:101];
%	end
%	
%	results.gridA  = -800000:5000:800000;
%	nA=numel(results.gridA);
%
%	results.thkA=NaN*ones(nA,length(results.timegrid));
%	results.sA=NaN*ones(nA,length(results.timegrid));
%	results.vxmeanA=NaN*ones(nA,length(results.timegrid));
%	results.vymeanA=NaN*ones(nA,length(results.timegrid));
%	results.maskA=NaN*ones(nA,length(results.timegrid));
%	
%	results.gridB  = -800000:5000:800000;
%	nB=numel(results.gridB);
%
%	results.thkB=NaN*ones(nB,length(results.timegrid));
%	results.sB=NaN*ones(nB,length(results.timegrid));
%	results.vxmeanB=NaN*ones(nB,length(results.timegrid));
%	results.vymeanB=NaN*ones(nB,length(results.timegrid));
%	results.maskB=NaN*ones(nB,length(results.timegrid));
%	
%	results.gridC  = -800000:5000:800000;
%	nC=numel(results.gridC);
%
%	results.thkC=NaN*ones(nC,length(results.timegrid));
%	results.sC=NaN*ones(nC,length(results.timegrid));
%	results.vxmeanC=NaN*ones(nC,length(results.timegrid));
%	results.vymeanC=NaN*ones(nC,length(results.timegrid));
%	results.maskC=NaN*ones(nC,length(results.timegrid));
%	
%	results.gridD  = -800000:5000:800000;
%	nD=numel(results.gridD);
%
%	results.thkD=NaN*ones(nD,length(results.timegrid));
%	results.sD=NaN*ones(nD,length(results.timegrid));
%	results.vxmeanD=NaN*ones(nD,length(results.timegrid));
%	results.vymeanD=NaN*ones(nD,length(results.timegrid));
%	results.maskD=NaN*ones(nD,length(results.timegrid));
%
%	results.gridE  = -800000:5000:800000;
%	nE=numel(results.gridE);
%
%	results.thkE=NaN*ones(nE,length(results.timegrid));
%	results.sE=NaN*ones(nE,length(results.timegrid));
%	results.vxmeanE=NaN*ones(nE,length(results.timegrid));
%	results.vymeanE=NaN*ones(nE,length(results.timegrid));
%	results.maskE=NaN*ones(nE,length(results.timegrid));
%
%	results.gridF  = -800000:5000:800000;
%	nF=numel(results.gridF);
%
%	results.thkF=NaN*ones(nF,length(results.timegrid));
%	results.sF=NaN*ones(nF,length(results.timegrid));
%	results.vxmeanF=NaN*ones(nF,length(results.timegrid));
%	results.vymeanF=NaN*ones(nF,length(results.timegrid));
%	results.maskF=NaN*ones(nF,length(results.timegrid));
%	
%	results.gridG  = -800000:5000:800000;
%	nG=numel(results.gridG);
%
%	results.thkG=NaN*ones(nG,length(results.timegrid));
%	results.sG=NaN*ones(nG,length(results.timegrid));
%	results.vxmeanG=NaN*ones(nG,length(results.timegrid));
%	results.vymeanG=NaN*ones(nG,length(results.timegrid));
%	results.maskG=NaN*ones(nG,length(results.timegrid));
%	
%	results.gridH  = -800000:5000:800000;
%	nH=numel(results.gridH);
%
%	results.thkH=NaN*ones(nH,length(results.timegrid));
%	results.sH=NaN*ones(nH,length(results.timegrid));
%	results.vxmeanH=NaN*ones(nH,length(results.timegrid));
%	results.vymeanH=NaN*ones(nH,length(results.timegrid));
%	results.maskH=NaN*ones(nH,length(results.timegrid));
%
%	%}}}
	
	%Create netcdf {{{
	disp(['creating netcdf...']);
	
	mode = netcdf.getConstant('NETCDF4');
	mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid=netcdf.create([directoryname '/CalvingMIP_' expnumber '_' modelname '_' institution '.nc'],mode);

	%General attributes
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Youngmin choi (yochoi@umd.edu)');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model','ISSM (Ice-sheet and Sea-level System Model)');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date',date());

	%Define dimensions
	ntime_id  = netcdf.defDim(ncid,'time',length(results.timescalar));
	nx_id  = netcdf.defDim(ncid,'x',length(results.gridx));
	ny_id  = netcdf.defDim(ncid,'y',length(results.gridy));
	
	disp(['creating netcdf (define variables)...']);
	%Define variables
	time_var_id = netcdf.defVar(ncid,'time','NC_FLOAT',[ntime_id]);
	
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
