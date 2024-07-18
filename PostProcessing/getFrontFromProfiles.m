function front=getFrontFromProfiles(p)
	% find exact ice front position
	x = p.distance;
	x0 = interpZeroPos(x, p.icemask)';
	% interpolate the other variables at x0
	fnames = fieldnames(p);
	for j = 1:numel(fnames)
		if (~strcmp(fnames{j}, 'distance'))
			temp = zeros(size(x0));
			for t = 1:numel(x0)
				temp(t) = interp1(x,p.(fnames{j})(:,t), x0(t));
			end
			front.(fnames{j}) = temp;
		end
	end
	front.distance = x0;
end
