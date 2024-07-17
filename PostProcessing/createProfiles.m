% Create csv files for the profiles used in CalvingMIP experiments
%		expnumber - experiment number [EXP1, ..., EXP5]
%		foldername - path to save the csv files
%		Nx - number of points, by default Nx = 200

function createProfiles(expnumber, foldername, Nx)

	flagCap = 0;
	flagHal = 0;
	flagCircle = 0;

	if (strcmp(expnumber,'EXP1') || strcmp(expnumber,'EXP2'))
		flagCircle = 1;
	elseif (strcmp(expnumber,'EXP3') || strcmp(expnumber,'EXP4'))
		flagCap = 1;
		flagHal = 1;
	elseif strcmp(expnumber,'EXP5')
		error([expnumber, ' not supported yet'])
	else
		error(['Unknown experiment ', expnumber])
	end

	if nargin < 3
		Nx = 200;
	end

	% set the profiles
	if flagCircle
		Circle.A.x = zeros(Nx, 1)*1e3;							Circle.A.y = linspace(0, 800, Nx)'*1e3;
		Circle.B.x = linspace(0, 800, Nx)'*1e3/sqrt(2);		Circle.B.y = Circle.B.x;
		Circle.C.x = linspace(0, 800, Nx)'*1e3;				Circle.C.y = zeros(Nx, 1)*1e3;
		Circle.D.x = linspace(0, 800, Nx)'*1e3/sqrt(2);		Circle.D.y = -Circle.D.x;

		Circle.E.x = zeros(Nx, 1)*1e3;							Circle.E.y = linspace(0, -800, Nx)'*1e3;
		Circle.F.x = linspace(0, -800, Nx)'*1e3/sqrt(2);	Circle.F.y = Circle.F.x;
		Circle.G.x = linspace(0, -800, Nx)'*1e3;				Circle.G.y = zeros(Nx, 1)*1e3;
		Circle.H.x = linspace(0, -800, Nx)'*1e3/sqrt(2);	Circle.H.y = -Circle.H.x;

		% prepare for table
		fn = fieldnames(Circle);
		Nfn = numel(fn);
		data = zeros(Nx, Nfn*3); % X, Y, S
		names = '';
		for i=1:Nfn
			data(:, i*3-2) = Circle.(fn{i}).x;
			names{i*3-2} = ['Circle_Profile_', upper(fn{i}), '_X'];
			data(:, i*3-1) = Circle.(fn{i}).y;
			names{i*3-1} = ['Circle_Profile_', upper(fn{i}), '_Y'];
			data(:, i*3) =  sqrt((Circle.(fn{i}).x-Circle.(fn{i}).x(1)).^2+(Circle.(fn{i}).y-Circle.(fn{i}).y(1)).^2);
			names{i*3} = ['Circle_Profile_', upper(fn{i}), '_S'];
		end
		% make a table
		T = array2table(data);
		T.Properties.VariableNames = names;
		writetable(T,[foldername, '/Circle_Profiles.csv'],'Delimiter',',','QuoteStrings',1)
	end
	% set the profiles
	if flagCap
		Caprona.A.x = linspace(-390, -590, Nx)'*1e3; Caprona.A.y = linspace(0, 450, Nx)'*1e3;
		Caprona.B.x = linspace(390, 590, Nx)'*1e3;   Caprona.B.y = linspace(0, 450, Nx)'*1e3;
		Caprona.C.x = linspace(-390, -590, Nx)'*1e3; Caprona.C.y = linspace(0, -450, Nx)'*1e3;
		Caprona.D.x = linspace(390, 590, Nx)'*1e3;   Caprona.D.y = linspace(0, -450, Nx)'*1e3;

		% prepare for table
		fn = fieldnames(Caprona);
		Nfn = numel(fn);
		data = zeros(Nx, Nfn*3); % X, Y, S
		names = '';
		for i=1:Nfn
			data(:, i*3-2) = Caprona.(fn{i}).x;
			names{i*3-2} = ['Caprona_Profile_', upper(fn{i}), '_X'];
			data(:, i*3-1) = Caprona.(fn{i}).y;
			names{i*3-1} = ['Caprona_Profile_', upper(fn{i}), '_Y'];
			data(:, i*3) =  sqrt((Caprona.(fn{i}).x-Caprona.(fn{i}).x(1)).^2+(Caprona.(fn{i}).y-Caprona.(fn{i}).y(1)).^2);
			names{i*3} = ['Caprona_Profile_', upper(fn{i}), '_S'];
		end
		% make a table
		T = array2table(data);
		T.Properties.VariableNames = names;
		writetable(T,[foldername, '/Caprona_Profiles.csv'],'Delimiter',',','QuoteStrings',1)
	end

	if flagHal
		Halbrane.A.x = ones(Nx,1)*(-150*1e3);        Halbrane.A.y = linspace(0, 740, Nx)'*1e3;
		Halbrane.B.x = ones(Nx,1)*(150*1e3);         Halbrane.B.y = linspace(0, 740, Nx)'*1e3;
		Halbrane.C.x = ones(Nx,1)*(-150*1e3);        Halbrane.C.y = linspace(0, -740, Nx)'*1e3;
		Halbrane.D.x = ones(Nx,1)*(150*1e3);         Halbrane.D.y = linspace(0, -740, Nx)'*1e3;

		% prepare for table
		fn = fieldnames(Halbrane);
		Nfn = numel(fn);
		data = zeros(Nx, Nfn*3); % X, Y, S
		names = '';
		for i=1:Nfn
			data(:, i*3-2) = Halbrane.(fn{i}).x;
			names{i*3-2} = ['Halbrane_Profile_', upper(fn{i}), '_X'];
			data(:, i*3-1) = Halbrane.(fn{i}).y;
			names{i*3-1} = ['Halbrane_Profile_', upper(fn{i}), '_Y'];
			data(:, i*3) =  sqrt((Halbrane.(fn{i}).x-Halbrane.(fn{i}).x(1)).^2+(Halbrane.(fn{i}).y-Halbrane.(fn{i}).y(1)).^2);
			names{i*3} = ['Halbrane_Profile_', upper(fn{i}), '_S'];
		end
		% make a table
		T = array2table(data);
		T.Properties.VariableNames = names;
		writetable(T,[foldername, '/Halbrane_Profiles.csv'],'Delimiter',',','QuoteStrings',1)
	end
