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

function J_history = costFunctionRate(X, y, theta, m, alpha, numberOfIterations)
  J_history = zeros(1, numberOfIterations);
  for i = 1:numberOfIterations      
    J_history(i) = costFunction(X, y, theta, m);     
    theta = getNewThetaValues(X, y, theta, m, alpha);       
  end
end

function theta = calculateTheta(X, y, theta, m, alpha)
  for i = 1:1500  
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

% Scale features and set them to zero mean
fprintf('Normalizing Features ...\n');

[X mu sigma] = featureNormalize(X);

######################################

x0 = ones(m, 1);
X = [x0, X ];

theta_init = zeros(1, columns(X));

fprintf('Calculating Cost Function... \n');
costFunction(X, y, theta_init, m);

alpha = 0.01;
theta = theta_init;

fprintf('Calculating theta... \n');
theta = calculateTheta(X, y, theta_init, m, alpha)
fprintf('Done \n');

figure;
plot(data(:, 1), y, 'rx', data(:, 1), X*theta', 'bx');
xlabel('Square ft');
ylabel('Price');

# Let's predict the price of a 1650 sq-ft, 3 br house
# X1 = 1650, X2 = 3
# At first I have to normalize these values, since theta was calculated 
# for normalized vals.
# X = (X .- mu) ./ sigma;

X = [1650, 3];
X = (X .- mu) ./ sigma;
X = [1, X];

fprintf('Predicted price is %f \n', X * theta');

hold on;
plot(1650, X * theta', 'k*');

# Now I wanna plot cost function changes for different learning rates
numberOfIterations = 50;
scale = 10^10;

figure;
plot(1:numberOfIterations, costFunctionRate(X, y, theta_init, m, 0.3, numberOfIterations) / scale);
xlabel('Number of iterations');
ylabel('Cost J (scaled)');
hold on;
plot(1:numberOfIterations, costFunctionRate(X, y, theta_init, m, 0.1, numberOfIterations) / scale);
hold on;
plot(1:numberOfIterations, costFunctionRate(X, y, theta_init, m, 0.03, numberOfIterations) / scale);
hold on;
plot(1:numberOfIterations, costFunctionRate(X, y, theta_init, m, 0.01, numberOfIterations) / scale);
legend ("alpha = 0.3", "alpha = 0.1", "alpha = 0.03", "alpha = 0.01");
