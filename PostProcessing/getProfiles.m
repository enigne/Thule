clear
close all

Cmap = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'};

% set the profiles
NC = 200;
Caprona.A.x = linspace(-390, -590, NC)'*1e3;	Caprona.A.y = linspace(0, 450, NC)'*1e3;
Caprona.B.x = linspace(390, 590, NC)'*1e3;	Caprona.B.y = linspace(0, 450, NC)'*1e3;
Caprona.C.x = linspace(-390, -590, NC)'*1e3; Caprona.C.y = linspace(0, -450, NC)'*1e3;
Caprona.D.x = linspace(390, 590, NC)'*1e3;	Caprona.D.y = linspace(0, -450, NC)'*1e3;
Halbrane.A.x = ones(NC,1)*(-150*1e3);			Halbrane.A.y = linspace(0, 740, NC)'*1e3;
Halbrane.B.x = ones(NC,1)*(150*1e3);			Halbrane.B.y = linspace(0, 740, NC)'*1e3;
Halbrane.C.x = ones(NC,1)*(-150*1e3);			Halbrane.C.y = linspace(0, -740, NC)'*1e3;
Halbrane.D.x = ones(NC,1)*(150*1e3);			Halbrane.D.y = linspace(0, -740, NC)'*1e3;

% project solutions to the profile
fn = fieldnames(Caprona);
for i=1:numel(fn)
	Caprona.(fn{i}) = project2Profile(md, Caprona.(fn{i}));
	Halbrane.(fn{i}) = project2Profile(md, Halbrane.(fn{i}));
end

% plot
figure
i = 1;
plot(Caprona.(fn{i}).xDist, Caprona.(fn{i}).bed, 'b')
hold on
for i=1:numel(fn)
	plot(Caprona.(fn{i}).xDist, Caprona.(fn{i}).base, 'color', Cmap{i})
	plot(Caprona.(fn{i}).xDist, Caprona.(fn{i}).base+Caprona.(fn{i}).thickness, 'color', Cmap{i})
end

figure
i = 1;
plot(Halbrane.(fn{i}).xDist, Halbrane.(fn{i}).bed, 'b')
hold on
for i=1:numel(fn)
	plot(Halbrane.(fn{i}).xDist, Halbrane.(fn{i}).base, 'color', Cmap{i})
	plot(Halbrane.(fn{i}).xDist, Halbrane.(fn{i}).base+Halbrane.(fn{i}).thickness, 'color', Cmap{i})
end

