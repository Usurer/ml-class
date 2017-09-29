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
  dJ = sum((h' - y) .* X) / m;

end

function thetas = getNewThetaValues(X, y, theta, m, alpha)
  delta = costFunctionDerivative(X, y, theta, m);
  thetas = theta - alpha * delta;
end

function theta = calculateTheta(X, y, theta, m, alpha)
  for i = 1:15
    theta = getNewThetaValues(X, y, theta, m, alpha);            
  end
end

function [X mu sigma] = featureNormalize(X)

  sigma = std(X);
  mu = sum(X) / length(X);
  X = (X .- mu) ./ sigma;  

end


data = load('ex1data2.txt');
X = data(:, 1:2);
y = data(:, 3);
m = length(y);

% Print out some data points
fprintf('First 10 examples from the dataset: \n');
fprintf(' x = [%.0f %.0f], y = %.0f \n', [X(1:10,:) y(1:10,:)]');

fprintf('Program paused. Press enter to continue.\n');
pause;

% Scale features and set them to zero mean
fprintf('Normalizing Features ...\n');

[X mu sigma] = featureNormalize(X)

fprintf('Normalizing done \n');

######################################

x1 = data(:,1);
y = data(:,2);

x0 = ones(m, 1);
theta_init = zeros(1,6);

X = [x0, x1, x1 .^ 2 / 10000, x1 .^ 3 / 10000, x1 .^ 4 / 10000, x1 .^ 5 / 10000 ];

costFunction(X, y, theta_init, m);

alpha = 0.01;
theta = theta_init;

theta = calculateTheta(X, y, theta_init, m, alpha)
plot(x1, y, 'rx', x1, theta * X', '-b');
#plot(x1, theta(1) + theta(2)*x1);

fprintf('Theta found by MY gradient descent: ');
fprintf('%f %f \n', theta(1), theta(2));

% Predict values for population sizes of 35,000 and 70,000
predict1 = theta(1) + 3.5*theta(2)
fprintf('For population = 35,000, we predict a profit of %f\n', predict1*10000);
predict2 = theta(1) + 7*theta(2)
fprintf('For population = 70,000, we predict a profit of %f\n', predict2*10000);
