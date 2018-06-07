function v = projball(x, y, r);

R = sqrt( sum(abs(x(:) - y(:)).^2) );
if R <= r
    v = x;
else
    v = y + (r/R) * (x-y);
end