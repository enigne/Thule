% 
% pf - struct of the profile, pf.x and pf.y are the coordinates

function pf = project2Profile(md, pf)
	px = pf.x;
	py = pf.y;
	pf.xDist = sqrt((px-px(1)).^2+(py-py(1)).^2);

	index = md.mesh.elements;
	x = md.mesh.x;
	y = md.mesh.y;
	pf.bed = InterpFromMeshToMesh2d(index, x, y, md.geometry.bed, px, py);
	pf.base = InterpFromMeshToMesh2d(index, x, y, [md.results.TransientSolution(:).Base], px, py);
	pf.thickness = InterpFromMeshToMesh2d(index, x, y, [md.results.TransientSolution(:).Thickness], px, py);
	pf.vx = InterpFromMeshToMesh2d(index, x, y, [md.results.TransientSolution(:).Vx], px, py);
	pf.vy = InterpFromMeshToMesh2d(index, x, y, [md.results.TransientSolution(:).Vy], px, py);
end
