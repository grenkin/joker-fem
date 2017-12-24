% MATLAB code

f = fopen('output.txt');
A = textscan(f, '%f %f %f %f');
x = A{1};
y = A{2};
phi = A{3};
theta = A{4};

K = 100;
N = K + 1;
X = zeros(N, N);
Y = X; Z = X; Zphi = X;
xi = 0;
i = 0;
for nX = 0 : K
    xi = xi + 1;
    yi = 0;
    for nY = 0 : K
        yi = yi + 1;
        i = i + 1;
        X(xi, yi) = x(i);
        Y(xi, yi) = y(i);
        Z(xi, yi) = theta(i);
        Zphi(xi, yi) = phi(i);
    end            
end

f2 = fopen('output_ff.txt');
A2 = textscan(f2, '%f %f %f %f');
x2 = A2{1};
y2 = A2{2};
phi2 = A2{3};
theta2 = A2{4};

K2 = 100;
N2 = K2 + 1;
X2 = zeros(N2, N2);
Y2 = X2; Z2 = X2;
xi = 0;
i = 0;
for nX = 0 : K2
    xi = xi + 1;
    yi = 0;
    for nY = 0 : K2
        yi = yi + 1;
        i = i + 1;
        X2(xi, yi) = x2(i);
        Y2(xi, yi) = y2(i);
        Z2(xi, yi) = theta2(i);
        Z2phi(xi, yi) = phi2(i);
    end            
end

figure
A = contour(X, Y, Z, 'k');
axis equal
clabel(A);

figure
A = contour(X, Y, Zphi, 'k');
axis equal
clabel(A);

%{
figure
A = contour(X2, Y2, Z2, 'k');
axis equal
clabel(A);
%}

D = abs(Z2 - Z);
Dphi = abs(Z2phi - Zphi);

max(max(D))
max(max(Dphi))