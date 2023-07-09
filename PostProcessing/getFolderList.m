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
		'20230705_EXP3_res_5000/',...
			};
		nameList = {
			'5000m, EXP3',...
			}; %}}}
	elseif Id == 4  % 4: EXP4 5km res{{{
		folderList = {
		'20230705_EXP3_res_5000/',...
		'20230705_EXP4_res_5000/',...
			};
		nameList = {
		'5000m, EXP 3',...
		'5000m, EXP 4',...
			}; %}}}
	end
