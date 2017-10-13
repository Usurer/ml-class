# Prevent Octave from thinking that this
# is a function file:

1;

function J = costFunction(X, y, theta, m)
  
  % I have issues with vector directions
  % So I have to transpose 'em from time to time
  % To get results I expect
  
  h = theta * X';
  J = sum((h' - y).^2) / (2*m);
  
end

function dJ = costFunctionDerivative(X, y, theta, m)

  h = theta * X';
  dJ = sum((h' - y) .* X)/m;

end

function thetas = getNewThetaValues(X, y, theta, m, alpha)
  delta = costFunctionDerivative(X, y, theta, m);
  l = 0;
  thetas = theta * (1 - alpha * l / m) - alpha * delta;
end

function theta = calculateTheta(X, y, theta, m, alpha)
  for i = 1:2000
    theta = getNewThetaValues(X, y, theta, m, alpha);
  end
end

function xNorm = normalize(X)

  xMin = min(X);
  xMax = max(X);
  xMean = sum(X) / length(X);

  xNorm = (X - xMean) / (xMax - xMin);

end

dataX = (load('ex5Linx.dat') + 1) ;
dataY = load('ex5Liny.dat');

plot(dataX, dataY, 'rx');

m = length(dataX);

xNormalized = normalize(dataX);

x0 = ones(m, 1);
x1 = xNormalized;
x2 = xNormalized .^2;
x3 = xNormalized .^3;
x4 = xNormalized .^4;
x5 = xNormalized .^5;
x6 = xNormalized .^6;
x7 = xNormalized .^7;

y = dataY;

theta_init = zeros(1,8);

X = [x0, x1, x2, x3, x4, x5, x6, x7]

fprintf('Cost function: ');
costFunction(X, y, theta_init, m)

fprintf('\n Cost function derivative: ');
costFunctionDerivative(X, y, theta_init, m)

alpha = 0.01;
theta = theta_init;


theta = calculateTheta(X, y, theta_init, m, alpha);

dataY;
theta * X';

plot(dataX, dataY, 'rx', dataX, theta * X', '-b');

