function b = calc_Kdotx(kcol, kval, diag,  x, b, nrrdof)

if 0
  for r=1:nrrdof
      b(r) = diag(r)*x(r);
  end

  for r=0:nrrdof-1
      lim = 10*r+kcol(10*r+1);
      for a = (10*r+2):lim
          b(r+1)= b(r+1) + kval(a)*x(kcol(a));
          b(kcol(a)) = b(kcol(a)) + kval(a)*x(r+1);
      end
  end
else
  b = diag .* x;
  for r = 1:nrrdof
    lim = kcol(r, 1);
    for a = 2:lim
      b(r) = b(r) + kval(r, a) * x(kcol(r, a));
      b(kcol(r, a)) = b(kcol(r, a)) + kval(r, a) * x(r);
    end
  end
end
end
