%% Time vs accuracy diagram for matrix exponential times vector computation

% close all
clearvars
% mydefaults

% Number of people in household - this also determines the matrix size
NN = [10]; % Each generates a matrix of size 66, 496 and 5050 respectively

% Loop for each household size inorder to vary the matrix size
for kk = 1:length(NN)
    
    % Number of people in a household
    N = NN(kk);
    
    % Set up transmission parameter within the household
    b = 0.1;
    alpha = 0.663;
    beta = b/((N-1)^alpha); %#ok<*PFOUS>
    gamma = 0.025;
    
    % Transmission between households
    % tau = 0.0047;
    tau = 0;
    
    % Create the generator matrix
    [Q,HHconfig] = SEI(N);
    
    % Generate the initial conditions vector
    tempI = find(HHconfig.dataI(:,3)==1); tempS = find(HHconfig.dataI(:,1)==N-1);
    pos = intersect(tempI,tempS);
    P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(pos,1) = 1;
    
    %
    Mfull = sparse(GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N));
    t = 360; % fixed time value - Interest is at t=360 days
    A = t*Mfull;
    exact = expm(full(A))*P0;  % TODO: do this in multiprecision (if A small)!
    lmin = full(min(diag(A)));
    lmax = full(max(diag(A)));
    runs = 100;
    
    %% DA1
    tim = []; err = [];
    %nsteps = 10:50:500;
    % Find the divisors of 360 (Time point of interest in days)
    nsteps = [1; unique(cumprod(perms(factor(360)),2))];
    for j = 1:length(nsteps),
        h = 360/nsteps(j);
        I = speye(size(A));
        B = I - h*Mfull;
        tic
        for k = 1:runs,
            appr = P0;
            for s = 1:nsteps(j),
                appr = B\appr;
            end
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    figure; set(gcf,'WindowStyle','Docked')
    loglog(1000*tim,err,'k-o')
    hold on
    
    %% DA2
    tim = []; err = [];
    %nsteps = 10:50:500;
    % Find the divisors of 360
    nsteps = [1; unique(cumprod(perms(factor(360)),2))];
    for j = 1:length(nsteps),
        h = 360/nsteps(j);
        I = speye(size(A));
        B = I - h*Mfull*(I - (h/2)*Mfull);
        tic
        for k = 1:runs,
            appr = P0;
            for s = 1:nsteps(j),
                appr = B\appr;
            end
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'k--s')
    
    %% DA3
    tim = []; err = [];
    %nsteps = 10:50:500;
    % Find the divisors of 360
    nsteps = [1; unique(cumprod(perms(factor(360)),2))];
    for j = 1:length(nsteps),
        h = 360/nsteps(j);
        I = speye(size(A));
        B = I - h*Mfull*(I - (h/2)*Mfull*(I - (h/3)*Mfull));
        tic
        for k = 1:runs,
            appr = P0;
            for s = 1:nsteps(j),
                appr = B\appr;
            end
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'k-.x')
    
    %% RK4
    tim = []; err = [];
    %nsteps = 10:50:500;
    % Find the divisors of 360
    nsteps = [1; unique(cumprod(perms(factor(360)),2))];
    for j = 1:length(nsteps),
        h = 360/nsteps(j);
        I = speye(size(A));
        B = I + h*Mfull*(I + (h/2)*Mfull*(I + (h/3)*Mfull*(I + (h/4)*Mfull)));
        tic
        for k = 1:runs,
            appr = P0;
            for s = 1:nsteps(j),
                appr = B*appr;
            end
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'r-x')
    % ylim([1e-16,1])
    
    %% EXPM
    tim = []; err = [];
    tic
    for k = 1:runs,
        appr = expm(A)*P0;
    end
    tim = toc/runs;
    err = norm(appr-exact,inf)/norm(exact,inf);
    semilogy(1000*tim,err,'b-+','MarkerSize',16)
    
    %% SEXPM
    tim = []; err = [];
    tic
    for k = 1:runs,
        appr = sexpmv(A,P0,0,'p');
    end
    tim = toc/runs;
    err = norm(appr-exact,inf)/norm(exact,inf);
    semilogy(1000*tim,err,'b-x','MarkerSize',16)
    
    %% EXPMV_2010
    tim = []; err = [];
    prec = {'half','single','double'};
    for j = 1:length(prec),
        tic
        for k = 1:runs,
            appr = expmv_2010(360,Mfull,P0,[],prec{j});
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'b:*','MarkerSize',16)
    
    %% EXPOKIT
    tim = []; err = [];
    tol = 10.^-(1:7);
    for j = 1:length(tol),
        tic
        for k = 1:runs,
            appr = mexpv(360, Mfull, P0, tol(j));
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'c--*')
    
    %% Chebyshev for exp(A)*P0 to various accuracies
    tim = []; err = [];
    tol = 10.^-(2:13);
    mmax = 50:10:80;
    for j = 1:length(mmax),
        tic
        for k = 1:runs,
            %[appr,st] = polycheby2(A, P0, tol(j), 4000, lmin, lmax); st
            [appr,st] = polycheby2(A, P0, 0, mmax(j), lmin, lmax);
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'g-*')
    
    %% rational Krylov (single shift)
    tim = []; err = [];
    mmax = 1:20; shift = 1;
    for j = 1:length(mmax), m = mmax(j);
        tic
        for k = 1:runs,
            appr = sikrylov(A,P0,m);
        end
        tim(j) = toc/runs;
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'m-+')
    
    %% ODE45 - Matlab's optimised one
    tim = []; err = [];
    tol = 10.^-(1:7);
    f = @(t,x)sparse(GenMatrixCalc(Q,beta,tau,gamma,HHconfig,x,N))*x;
    for j = 1:length(tol),
        tic
        for k = 1:runs,
            [t1,appr1] = ode45(f,[0 360],P0,odeset('RelTol',tol(j)));%
        end
        tim(j) = toc/runs;
        appr = appr1(end,:); appr = appr';
        err(j) = norm(appr-exact,inf)/norm(exact,inf);
    end
    loglog(1000*tim,err,'g-+')
    
    %% Some plotting to do
    ylim = get(gca,'YLim'); set(gca,'YLim',[1e-16 10^11],'XLim',[10^-3 10^6])
    xlabel('computational time (ms)'); ylabel('accuracy')
    hand = legend('DA1','DA2','DA3','RK4','EXPM','SEXPM','EXPMV2010','EXPOKIT','CHEBY','RATKRYL','ODE45');
    hand.FontSize = 12;
    set(hand,'Box','off')
    mypdf('workacc',.8,2)
    
end