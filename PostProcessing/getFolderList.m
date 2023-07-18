function [folderList, nameList] = getFolderList(Id)
	%getFolderList - return the working folders
	%
	% Author: Cheng Gong
	% Last modified: 2023-06-19
	if nargin < 1
		Id = 0;
	end
	if Id == 0  % Latest test{{{
		folderList = {
		'20230714_EXP4_res_2500/',...
			};
		nameList = {
			'2500m, EXP4',...
			}; %}}}
	elseif Id == 4  % 4: EXP3&4 5km res{{{
		folderList = {
		'20230705_EXP3_res_5000/',...
		'20230709_EXP4_res_5000/',...
			};
		nameList = {
		'5000m, EXP 3',...
		'5000m, EXP 4',...
			}; %}}}
	elseif Id == 5  % 5: EXP3&4 2.5km res{{{
		folderList = {
		'20230714_EXP3_res_2500/',...
		'20230714_EXP4_res_2500/',...
			};
		nameList = {
		'2500m, EXP 3',...
		'2500m, EXP 4',...
			}; %}}}
	end
