function[ x objval x2 time ] = Lshap( s , maxrun )
% Solving two stage stochastic programing using Lshape algorithm
%--------------------Input example-----------------------------
% % ----------------- stage 1 -----------------
% % coef in the obj function for decision variable in stage 1
% smodel.x1.coef = [100 150];
% smodel.x1.start = [];
% smodel.x1.lb = [40 20];
% smodel.x1.ub = [inf inf];
% smodel.x1.ctype ='CC';
% 
% %constrains for x1 only
% smodel.x1.lhs = -inf;
% smodel.x1.cos = [ 1 1 ];
% smodel.x1.rhs = 120;
% 
% 
% 
% 
% % the coef in the constrains, noticing that in the stochastic programing,
% 
% % ----------------- stage 2 -----------------
% % determined coef in the obj function for decision variable in stage 2
% smodel.x2.lb = [0  0];
% smodel.x2.ub = [inf inf];
% smodel.x2.ctype ='CC';
% 
% % probabilities corresponding to each scenario
% smodel.x2.prob = [0.4 0.6];
% % parameters under each scenario
% % [q1 q2]
% smodel.x2.coef = [-24 -28 ;-28 -32];
% % [d1 d2]
% smodel.x2.h = [500 100; 300 300];
% 
% %constrains for x2
% smodel.x2.cos = [ 6 10; 8 5; 1 0; 0 1];
% 
% % rhs
% smodel.x2.rhs.c = [0 0 0 0]';
% smodel.x2.rhs.d = [0 0; 0 0; 1 0; 0 1];
% % x1 are on the rh side
% % noticing that the coef of x1 are static here in stage 2
% smodel.x2.rhs.x1 = [60 0; 0 80 ; 0 0; 0 0];
% [ rslt objval x2 ] = Lshap(smodel);
%---------------------------------------------------------------
    tic;
    if nargin==1;
        maxrun = 100;
    end
    x1_ = s.x1.start;
    
    % in case of accuration problem
    IPSILON = 0.0000001;
    
    mStage1 = Cplex('stage1');
    mStage1.Model.sense = 'minimize';
    mStage1.addCols(s.x1.coef',[],s.x1.lb',s.x1.ub',s.x1.ctype);
    mStage1.addRows(s.x1.lhs,s.x1.cos,s.x1.rhs);

    if size(x1_) == [0 0]
        fprintf('start point calculating\n');
        % calculation of the start x1_    
        % Solve model
        mStage1.solve();
        x1_=mStage1.Solution.x;
    else
        fprintf('start point given\n');
    end
    
    % the number of desicion variable in 1st stage
    x1n = size(s.x1.lb,2);
    
    % add theta
    mStage1.addCols( 1 , [] ,  -inf, inf );

    % the number of desicion variable in 2nd stage
    x2n = size(s.x2.lb,2);
    % the number of scenario
    scnrn = size(s.x2.prob);
    scnrn = scnrn(2);

    models = cell(1,scnrn);
    for  sc = 1 : scnrn
        m = Cplex('model');
        m.Model.sense = 'minimize';
        m.addCols(s.x2.coef(sc,:)',[],s.x2.lb',s.x2.ub');
        % TODO:Integer!!!!!!!!!!!!!!!!!!!
        % m.addCols(s.x2.coef(sc,:)',[],s.x2.lb',s.x2.ub',s.x2.ctype);
        models{sc} = m;
    end


    % Loop start
    theta = -inf;
    E = zeros(scnrn,1);
    e = 0;
    flag = 0;
    x2_ = zeros( scnrn, x2n );
    
    for i = 1 : maxrun
        fprintf('\n----------------------Iteration %d:----------------------\n',i)    
        if theta > -inf 
            % add a cut and calculate new x1_
            mStage1.addRows(e, [E 1] ,inf);
            mStage1.solve();
            x1_=mStage1.Solution.x(1:end-1);
            theta = mStage1.Solution.x(end);
        end
        fprintf('\nx1_:')
        disp (x1_');
        esum = zeros(scnrn,1);
        Esum = zeros(scnrn,x1n);
        for sc = 1 : scnrn
            fprintf('\n Scenario %d:\n',sc)
            m = models{sc};
            consn = size(m.Model.A);
            consn = consn(1);
            if consn >0
                %remove the constrains
                m.delRows([1:consn]);
            end
            rhs = s.x2.rhs.x1{sc}*x1_ +  s.x2.rhs.h{sc};
            lhs = rhs;
            lhs( s.x2.rhs.equation == 0 ) = -Inf;
            cos = s.x2.cos;
            m.addRows(lhs,cos,rhs);
            m.solve();
            fprintf('\nObjective value: %d\n',m.Solution.objval);
            h = s.x2.rhs.h{sc};
            esum(sc,:) = m.Solution.dual'*h;
            T = - s.x2.rhs.x1{sc};
            Esum(sc,:) = m.Solution.dual'*T;
            x2_(sc,:) = m.Solution.x;
        end
        % calculating theta
        e = s.x2.prob * esum;
        E = s.x2.prob * Esum;
        w = e - E*x1_;
        fprintf('\n-----------------#End of Iteration %d:---------------\n',i)
    
        if theta < w && w - theta > IPSILON 
            fprintf('add a cut theta=%.2f; e=%.2f; E=',w,e);
            disp(E);
            fprintf('\n');
            theta = w;
        else
            fprintf('Optimal solution found in Iteration %d\n',i);
            flag = 1;
            break;
        end
    end
    if flag == 0 
        fprintf('Number of iteration exceeded the limit\n');
    end
    
    %return fields
    x = x1_;
    objval = mStage1.Solution.objval;
    x2 = x2_;
    toc;
    time = toc;
end