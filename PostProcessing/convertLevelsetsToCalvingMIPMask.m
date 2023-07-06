function mask = convertLevelsetsToCalvingMIPMask(ice_levelset, ocean_levelset)

   mask = ice_levelset;
   mask(ice_levelset<0) = 1;     % grounded and floating
   mask(ice_levelset>0) = 3;     % open ocean
   mask((ocean_levelset<0) & (ice_levelset<0)) = 2;

