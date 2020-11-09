function [fig_return] = show_rosenbrock()
nnx = 100;
x1 = linspace(-2.048, 2.048, nnx);
x2 = x1;
[X1, X2] = meshgrid(x1,x2);
Z = 100 .* (X2 - X1.^2).^2 + (1 - X1).^2;
fig_return = meshc(X1,X2, Z);
colorbar;
view([-211 20]);
end