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


dataX = load('ex5Linx.dat');
dataY = load('ex5Liny.dat');

plot(dataX, dataY, 'rx');

m = length(dataX);

x0 = ones(m, 1);
x1 = dataX;
x2 = dataX .^2;
x3 = dataX .^3;
x4 = dataX .^4;
x5 = dataX .^5;
x6 = dataX .^6;
x7 = dataX .^7;

y = dataY;

theta_init = zeros(1,8);

X = [x0, x1, x2, x3, x4, x5, x6, x7]

fprintf('Cost function: ');
costFunction(X, y, theta_init, m)

fprintf('\n Cost function derivative: ');
costFunctionDerivative(X, y, theta_init, m)

alpha = 0.01;
theta = theta_init;


theta = calculateTheta(X, y, theta_init, m, alpha)

dataY
theta * X'

plot(dataX, dataY, 'rx', dataX, theta * X', '-b');

