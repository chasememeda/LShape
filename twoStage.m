%%Two stage

smodel.type = 'two stage';
% ----------------- stage 1 -----------------
% coef in the obj function for decision variable in stage 1
smodel.x1.coef = [100 150];
smodel.x1.start = [];
smodel.x1.lb = [40 20];
smodel.x1.ub = [inf inf];
smodel.x1.ctype ='CC';

%constrains for x1 only
smodel.x1.lhs = -inf;
smodel.x1.cos = [ 1 1 ];
smodel.x1.rhs = 120;




% the coef in the constrains, noticing that in the stochastic programing,

% ----------------- stage 2 -----------------
% determined coef in the obj function for decision variable in stage 2
smodel.x2.lb = [0  0];
smodel.x2.ub = [inf inf];
smodel.x2.ctype ='CC';

% probabilities corresponding to each scenario
smodel.x2.prob = [0.4 0.6];
% parameters under each scenario
% [q1 q2]
smodel.x2.coef = [-24 -28 ;-28 -32];

%constrains for x2
smodel.x2.cos = [ 6 10; 8 5; 1 0; 0 1];

% rhs
smodel.x2.rhs.h = { [0,0,500,100]';[0,0,300,300]'};
smodel.x2.rhs.equation = [ 0 0 0 0 ]';
% x1 are on the rh side
% noticing that the coef of x1 are static here in stage 2
%smodel.x2.rhs.x1 = {[60 0; 0 80 ; 0 0; 0 0];[60 0; 0 80 ; 0 0; 0 0]};
smodel.x2.rhs.x1 = {sparse([60 0; 0 80 ; 0 0; 0 0]);sparse([60 0; 0 80 ; 0 0; 0 0])};

[ rslt objval x2 ] = Lshap(smodel);
