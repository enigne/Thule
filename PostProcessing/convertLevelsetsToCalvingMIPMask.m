function mask = convertLevelsetsToCalvingMIPMask(ice_levelset, ocean_levelset)
	%convertLevelsetsToCalvingMIPMask:
	%	To convert levelset functions defined ice domain to the mask for CalvingMIP
	%	where grounded=1, floating=2, open ocean=3, 
	%
	%	Usage:
	%		mask = convertLevelsetsToCalvingMIPMask(ice_levelset, ocean_levelset);
	%
	%  ice_levelset and ocean_levelset are the levelset functions for ice and ocean as defined in ISSM
	%
	%	Example:
	%		mask = convertLevelsetsToCalvingMIPMask(md.mask.ice_levelset, md.mask.ocean_levelset);
	%		mask = convertLevelsetsToCalvingMIPMask(md.results.TransientSolution(end).MaskIceLevelset, md.results.TransientSolution(:).MaskOceanLevelset);

	mask = ice_levelset;
	mask(ice_levelset <= 0) = 1;    % ice grounded and floating
	mask(ice_levelset > 0) = 3;     % open ocean / no ice
	mask((ocean_levelset <= 0) & (ice_levelset <= 0)) = 2;

