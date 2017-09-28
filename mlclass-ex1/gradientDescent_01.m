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

function delta = getDelta(X, y, theta, m, j)

  h = theta * X';
  delta = sum((h' - y)' * X(:,j)) / m;

end

function theta_j = getThetaJ(X, y, theta, m, j, alpha)
  delta = getDelta(X, y, theta, m, j);
  theta_j = theta(j) - alpha * delta;
end

# I want to do the same but for all thetas in one function, using matrix

# At first I'll get vector of dJ/dTheta, i.e. cost function derivatives
# My fist version (see getDelta) didn't use matrix properly and index of X was required
function deltas = getDeltas(X, y, theta, m)

  h = theta * X';
  deltas = sum((h' - y) .* X) / m;

end

function thetas = getThetas(X, y, theta, m, alpha)
  delta = getDeltas(X, y, theta, m);
  thetas = theta - alpha * delta;
end

function theta = getThetas2(X, y, theta, m, alpha)
  for i = 1:1500
    theta = getThetas(X, y, theta, m, alpha);        
    #costFunction(X, y, theta_new, m);
  end
end


# End new matrix thetas calculation

function theta = getTheta(X, y, theta, m, alpha)
  for i = 1:1500
    theta_1 = getThetaJ(X, y, theta, m, 1, alpha);
    theta_2 = getThetaJ(X, y, theta, m, 2, alpha);
    theta(1) = theta_1;
    theta(2) = theta_2;
    costFunction(X, y, theta, m);
  end
end

data = load('ex1data1.txt');
m = length(data);

x1 = data(:,1);
y = data(:,2);

x0 = ones(m, 1);
theta_init = zeros(1,2);

X = [x0, x1];

costFunction(X, y, theta_init, m);

alpha = 0.01;
theta = theta_init;

theta = getTheta(X, y, theta_init, m, alpha)
thetas = getThetas2(X, y, theta_init, m, alpha)
%plot(x1, y, 'rx', x1, theta(1) + theta(2)*x1);
%plot(x1, theta(1) + theta(2)*x1);

fprintf('Theta found by MY gradient descent: ');
fprintf('%f %f \n', theta(1), theta(2));

% Predict values for population sizes of 35,000 and 70,000
predict1 = theta(1) + 3.5*theta(2)
fprintf('For population = 35,000, we predict a profit of %f\n', predict1*10000);
predict2 = theta(1) + 7*theta(2)
fprintf('For population = 70,000, we predict a profit of %f\n', predict2*10000);
