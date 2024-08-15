function B = set_matrix_B(B,r,s,def)
	VOXSIZE = def.VOXSIZE;
	% r,s are the local coordinates in the isoparametric element
	% constants in shape functions for node 0 to 7
	 kr = [-1,  1,  1, -1];
	 ks= [-1, -1,  1,  1];
	% shape function and derivatives in point r,s

	for (i=1:4)   % for all nodes

		% value of the shape function belonging to node i
		% N[i] = .25 * (1 + kx[i] * localx) * (1 + ky[i] * localy);
		% dNdr     dxdr dydr      dNdx      dxdr   0        dNdx
		% dNds  =  dxds dyds   *  dNdy  =     0  dyds    *  dNdy

		% rewriting gives the values of the derivatives of the shape function
		dNdx(i) = (2/VOXSIZE) * .25 * kr(i) * (1 + ks(i) * s);
		dNdy(i) = (2/VOXSIZE) * .25 * ks(i) * (1 + kr(i) * r);

		% calculate strain displacement matrix B
		B(1, (2*i-1)) 	 = dNdx(i);
		B(1, (2*i)) =       0;
		B(2,  (2*i-1)) 	 =       0;
		B(2, (2*i)) = dNdy(i);
		B(3,  (2*i-1))	 = dNdy(i);
		B(3, (2*i)) = dNdx(i);
	end % endfor all nodes
end
