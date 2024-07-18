function profiles = project2Profiles(md, P, suffixname)
	nameList = 'A':char('A'+size(P,2)/3-1);
	% get the model solution into a matrix
	x = md.mesh.x;
	y = md.mesh.y;
	index = md.mesh.elements;

	% extract transient solutions
	profiles.time = [md.results.TransientSolution(:).time];
	vx = [md.results.TransientSolution(:).Vx];
	vy = [md.results.TransientSolution(:).Vy];
	thickness = [md.results.TransientSolution(:).Thickness];
	icemask = [md.results.TransientSolution(:).MaskIceLevelset];
	oceanmask = [md.results.TransientSolution(:).MaskOceanLevelset];
	

	for i = 1:numel(nameList)
		disp(['  Projecting solutions onto ' suffixname, ' profile ', nameList(i)])
		pfx = P.([suffixname, '_Profile_', nameList(i), '_X']);
		pfy = P.([suffixname, '_Profile_', nameList(i), '_Y']);

		% use a temporary profile variable to hold all the info
		pf = [];
		% distance from the start
		pf.distance = P.([suffixname, '_Profile_', nameList(i), '_S']);

		% project ice thickness, icemask, vx and vy
		pf.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, pfx, pfy);
		pf.vx = InterpFromMeshToMesh2d(index, x, y, vx, pfx, pfy);
		pf.vy = InterpFromMeshToMesh2d(index, x, y, vy, pfx, pfy);
		pf.icemask = InterpFromMeshToMesh2d(index, x, y, icemask, pfx, pfy);
		pf.oceanmask = InterpFromMeshToMesh2d(index, x, y, oceanmask, pfx, pfy);
		pf.mask = convertLevelsetsToCalvingMIPMask(pf.icemask, pf.oceanmask);
		pf.vel = sqrt(pf.vx.^2+pf.vy.^2);

		profiles.([suffixname, '_', nameList(i)]) = pf;
	end% pf - struct of the profile
end
