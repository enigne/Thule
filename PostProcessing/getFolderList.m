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
		'20230616_EXP3_res_2000_relaxT_1000/',...
			'20230616_EXP3_res_5000/',...
			};
		nameList = {
		'2000m, relax 1000',...
			'5000m, relax 1000',...
			}; %}}}
	elseif Id == 4  % 4: EXP4 5km res{{{
		folderList = {
		'20230623_EXP4_res_5000/',...
			};
		nameList = {
		'5000m, EXP 4',...
			}; %}}}
	end
