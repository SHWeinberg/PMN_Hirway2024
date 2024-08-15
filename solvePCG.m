function u = solvePCG(kcol, kval, u, f, nrrdof, def)
ACCURACY = def.ACCURACY;
if 0
  ui = zeros(nrrdof,1);
  for i=1:nrrdof
      ui(i)=u(i);
  end
  ri   = zeros(nrrdof,1);
  diag = zeros(nrrdof,1);
  invC = zeros(nrrdof,1);
  zi   = zeros(nrrdof,1);
  pi   = zeros(nrrdof,1);
  qi   = zeros(nrrdof,1);

  for i=0:nrrdof-1
      diag(i+1) = kval(10*i+1);
  end

  for i=1:nrrdof % for each row in K
      if diag(i) ~= 0.0 % if Kii not zero
          invC(i) = 1.0/diag(i);  					% invC = inv(diag(K))
      else
          invC(i) = 0.0;
      end
  end

  qi = calc_Kdotx(kcol, kval, diag, ui, qi, nrrdof);

  for i=1:nrrdof
      ri(i) = f(i)-qi(i); 								% r0 = f-K*u0
  end

  for i=1:nrrdof
      zi(i) = invC(i)*ri(i); 							% z0 = inv(C)*r0
  end

  for i=1:nrrdof
      pi(i) = zi(i); 									% p0 = z0
  end

  rhoinew = .0;
  for i=1:nrrdof
      rhoinew = rhoinew + ri(i)*zi(i); 						% rhoi = zi*ri
  end
  initrho = .0;
  for i=1:nrrdof
      initrho = initrho + invC(i)*f(i)*f(i);					% FOR ACCURACY
  end
  % start iterative solve
  iter = 0;
  % delta(iter) = rhoinew-ACCURACY*initrho;
  while rhoinew>ACCURACY*initrho
      % progress_total = progress_total + 0.05 + 0.6 + 0.25*(iter);
      % progresbar(progress_total,[],[],iter,[])
      rhoi = rhoinew;

      qi = calc_Kdotx(kcol,kval,diag, pi, qi, nrrdof); 	% qi = K*pi

      pq = 0;
      for i=1:nrrdof
          pq = pq + pi(i)*qi(i);
      end

      alfi = rhoi/pq;									% alfi = rhoi/(pi*qi)

      for i=1:nrrdof
          ui(i)= ui(i) + alfi*pi(i);							% ui+1 = ui+alfi*pi
      end

      for i=1:nrrdof
          ri(i) = ri(i) - alfi*qi(i);							% ri+1 = ri-alfi*qi
      end

      for i=1:nrrdof
          zi(i) = invC(i) * ri(i); 						% zi+1 = inv(C)*ri+1
      end

      rhoinew = .0;
      for i=1:nrrdof
          rhoinew = rhoinew + ri(i) * zi(i); 					% rhoi+1 = ri+1*zi+1
      end

      beti = rhoinew / rhoi;							% beti = rhoinew/rhoi

      for i=1:nrrdof
          pi(i) = zi(i) + beti * pi(i); 				% pi+1 = zi+1 + betai*pi
      end

      % if mod(iter,10) == 0
      %     fprintf('i %4d, rhoinew/initrho=%18.11d\n',iter, rhoinew/initrho);
      % end
      iter = iter + 1;
  end
  % fprintf('Stop iterating at iter %d\n',iter);
  for i=1:nrrdof
      u(i) = ui(i);
  end
else
  ui = u;
  ri(nrrdof, 1) = 0;
  diag = ri; invC = ri; zi = ri; pi = ri; qi = ri;

  diag = kval(:, 1);
  for i = 1:nrrdof % for each row in K
      if diag(i) ~= 0.0 % if Kii not zero
          invC(i, 1) = 1.0/diag(i);  					% invC = inv(diag(K))
      else
          invC(i, 1) = 0.0;
      end
  end
  % invC(diag ~= 0) = 1/diag(diag ~= 0);

  qi = calc_Kdotx(kcol, kval, diag, ui, qi, nrrdof);

  ri = f - qi;
  zi = invC .* ri;
  pi = zi;

  rhoinew = .0;
  for i=1:nrrdof
      rhoinew = rhoinew + ri(i)*zi(i); 						% rhoi = zi*ri
  end
  initrho = .0;
  for i=1:nrrdof
      initrho = initrho + invC(i) * f(i) * f(i);					% FOR ACCURACY
  end
  % start iterative solve
  iter = 0;
  % delta(iter) = rhoinew-ACCURACY*initrho;
  while rhoinew>ACCURACY*initrho
      % progress_total = progress_total + 0.05 + 0.6 + 0.25*(iter);
      % progresbar(progress_total,[],[],iter,[])
      rhoi = rhoinew;

      qi = calc_Kdotx(kcol,kval,diag, pi, qi, nrrdof); 	% qi = K*pi

      pq = 0;
      for i = 1:nrrdof
          pq = pq + pi(i) * qi(i);
      end

      alfi = rhoi/pq;									% alfi = rhoi/(pi*qi)

      for i = 1:nrrdof
          ui(i)= ui(i) + alfi*pi(i);							% ui+1 = ui+alfi*pi
      end

      for i=1:nrrdof
          ri(i) = ri(i) - alfi * qi(i);							% ri+1 = ri-alfi*qi
      end

      zi = invC .* ri;

      rhoinew = .0;
      for i = 1:nrrdof
          rhoinew = rhoinew + ri(i) * zi(i); 					% rhoi+1 = ri+1*zi+1
      end

      beti = rhoinew / rhoi;							% beti = rhoinew/rhoi

      for i = 1:nrrdof
          pi(i) = zi(i) + beti * pi(i); 				% pi+1 = zi+1 + betai*pi
      end

      % if mod(iter,10) == 0
      %     fprintf('i %4d, rhoinew/initrho=%18.11d\n',iter, rhoinew/initrho);
      % end
      iter = iter + 1;
  end

  % fprintf('Stop iterating at iter %d\n',iter);
  u = ui;

end
end
