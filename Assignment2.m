%% ELEC 4700 Assignment 2: Finite Difference Method
% Student Name: Cyril Onyebuchie
% (ID: 101133448)

%% Introduction
% The MATLAB code is design to perform the Finite Difference Method to solve
% electrostatic potential in rectangular shape region with dimensions L*W 
% where dimensions ratio (L/W) is 3/2. 
% Firstly, the case V=Vo at x=0 & V=0 at x=L. The case is solved using
% MATLAB code as given below. 
% Then, the next parts are solved for other parameter such as V=Vo at x=0,
% x=L and V=0 at y=0, and y=W.
% All these results are calculated with analytical solution. The matrix
% form is used to solve the problems. The variables are assigned for the 
% boundary condtions so that it can be changed without any major change in
% the code. 

%% Question 1

clear all
clearvars -GLOBAL
format shortE
close all

% Input Parameters
sz = 20; % Rectangular size
Vbo = 1; % Voltage boundary conditions
Ln = 3*sz; % Recangle length
Wd = 2*sz; % Rectangle width
GM1=sparse(1, Ln); % Matrix to solve problem through iterations method
GM2=sparse(Wd, Ln); % Matrix to solve problem through iterations method
GM3=sparse(Wd, Ln); % Matrix to solve problem through Analytical method
BoCo1 = 0;      % Q1 Boundary condition
BoCo2 = Vbo;    % Q1 Boundary condition

iter = 100; % Number of iterations

%% Question 1 - Part (a)

fig1 = figure; 
GM1(Ln)= BoCo1; % Seting of boundary conditions (V=0 at x=L)
GM1(1)= BoCo2;  % Setting of boundary conditions (V=Vo at x=0)
GM1(2:Ln-1) = deal(0.5*(BoCo1+BoCo2));

figure(fig1)   % Plot figure

for steps =1:iter
    plot(GM1);
    title('1-D finite Difference Process')
    pause(0.01);
    for ite=2:Ln-1
        GM1(ite)= 1/2*(GM1(ite-1)+GM1(ite+1));
    end
end

%% Question 1 - Part (b)
% Now, the part of question 1 has been solved. Firstly, the analytical
% solution is done and then compared with the meshing. 

Ntrm=iter;
fig2=figure;

figure(fig2)
for ite = 1:Ln       
    for jd = 1:Wd
        GM3(jd,ite) = AnalyticalSolution(ite,jd,Vbo,Ln,Wd,Ntrm); 
    end
end
subplot(3,1,1);
surf(GM3);
title('Analytical Solution');
pause(0.01);

% Now, the 2D matrix is populated for new boundary conditions
% V = Vo at x = 0, x = L and V = 0 at y = 0, y = W

% For the average boundry conditions, populated
[GM2(1, 1), GM2(1, Ln), GM2(Wd, 1), GM2(Wd, Ln), ...  
    GM2(2 : Wd - 1, 2 : Ln - 1)] = deal (1/2*(BoCo1 + BoCo2));
[GM2(1, 2 : Ln-1), GM2(Wd, 2 : Ln-1)] = deal(BoCo1);
[GM2(2 : Wd-1, 1), GM2(2 : Wd-1, Ln)] = deal(BoCo2);

for step = 1:iter
    subplot(3,1,2);
    surf(GM2);
    title('Numerical Solution');
    subplot(3, 1, 3);
    surf(log(abs(GM2 - GM3)))
    set(gca, 'ZScale' , 'log')
    title('Difference in Numerical and Analytical Solution')
    pause(0.01);
    
    for ite = 2 : Ln-1
        for jd = 2:Wd-1
            GM2(jd, ite) = 0.25*(GM2(jd + 1, ite) + GM2(jd - 1, ite) + ...
                GM2(jd, ite + 1) + GM2(jd, ite - 1));
        end
    end
end

%% Conclusion Question 1
% It can be seen from the results that, the iteration method gives the 
% same result as the analytical method in 100 iterations. The error between
% the both methods is around -10E-2 that is error range to stop the
% analytical method. The results have also been tested for the iterations
% more than 100 but the results remain almost same without any big change.
% Thus, the numerical value is easier because this approximates to the
% actual result for the complex analytical problems. The process of meshing
% is also simple but it takes more time. The disadvantage of numerical
% method includes, this process is difficult for the non-rectangular and
% non-uniform geometric shapes. 

%% Question 2
% In this question, Finite Difference Method is used to solve current flow
% in rectangular shape having different resistivity and having higher
% resistivity outside as compared to the inside resitivity. The plots have
% been done for sigma(x, y), V(x, y) Ex, Ey, & J(x,y).

%% Input Parameters
% Bottle neck length is about 1/3 of total rectangular region length while 
% width is around 2/8 of total width.

cMp2 = ones(Ln,Wd); % Mapping for conductivity
LGp = Ln/3;         % Bottle nech length
WGp = 3*Wd/8;       % Width of boxes
sg = 10e-2;
[cMp2(floor((Ln - LGp)/2) : floor((Ln + LGp)/2), 1 : floor((Wd - WGp)/2)),...
    cMp2(floor((Ln - LGp)/2) : floor((Ln + LGp)/2), ...
    floor((Wd + WGp)/2) : Wd)] = deal(sg);

GM4 = sparse(Ln*Wd, Ln*Wd); 
BM3 = zeros(1, Ln*Wd);

[BM3(2 : Wd - 1), BM3((Ln - 1)*Wd + 1 : Ln*Wd- 1 )] = deal(BoCo2);
[BM3(1), BM3(Wd), BM3((Ln - 1)*Wd), BM3(Ln*Wd)] = deal(0.5*(BoCo1 + BoCo2)); 

for ite = 1 : Ln
    for jd = 1 : Wd
        n = jd + (ite - 1)*Wd;   
        if ite==1 || ite == Ln || jd == 1 || jd == Wd
            GM4(n,n)=1;
        else
            nxpp = jd + ite*Wd;
            nxmm = jd + (ite - 2)*Wd;
            nyp=n+1;
            nym=n-1;
            
            rxpp = 1/2*(cMp2(ite, jd) + cMp2(ite + 1, jd));
            rxmm = 1/2*(cMp2(ite, jd) + cMp2(ite - 1, jd));
            rypp = 1/2*(cMp2(ite, jd) + cMp2(ite, jd + 1));
            rymm = 1/2*(cMp2(ite, jd) + cMp2(ite, jd - 1));
            
            GM4(n, n) = -(rxpp + rxmm + rypp + rymm);
            GM4(n, nxpp) = rxpp;
            GM4(n, nxmm) = rxmm;
            GM4(n, nyp) = rypp;
            GM4(n, nym) = rymm;
        end    
    end
end

V2V=GM4\(BM3'); % Variable to get V(x,y) from conductivity map

V2Mp = zeros(Ln, Wd);

for ite = 1 : Ln
    for jd = 1 : Wd
        V2Mp(ite, jd) = V2V(jd + (ite - 1)*Wd);
    end
end

eXs = zeros(Ln, Wd);
eYs = zeros(Ln, Wd);

for ite = 1 : Ln
    for jd = 1 : Wd
        if ite == 1
            eXs(ite,jd) = V2Mp(ite + 1, jd) - V2Mp(ite, jd);
        elseif ite == Ln
            eXs(ite, jd) = V2Mp(ite,jd) - V2Mp(ite - 1, jd);
        else
            eXs(ite, jd) = (V2Mp(ite + 1, jd) - V2Mp(ite - 1, jd))*0.5;
        end
        if jd == 1
            eYs(ite,jd) = V2Mp(ite,jd+1)-V2Mp(ite,jd);
        elseif jd == Wd
            eYs(ite,jd) = V2Mp(ite,jd)-V2Mp(ite,jd-1);
        else
            eYs(ite,jd) = (V2Mp(ite,jd+1)-V2Mp(ite,jd-1))*0.5;
        end
    end
end

eXs=-eXs;
eYs=-eYs;

%% Current Density Vectors
% The current density vecots are designed.
JxD = cMp2.*eXs; % X value for Current Density
JyD = cMp2.*eYs; % Y value for Current Density 

pointSampl = 3;
[x, y]= meshgrid(1:floor(Wd/pointSampl),1:floor(Ln/pointSampl));
x = pointSampl*x - pointSampl + 1;
y = pointSampl*y - pointSampl + 1;

eXSamp = zeros(floor(Ln/pointSampl),floor(Wd/pointSampl));
eYSamp = zeros(floor(Ln/pointSampl),floor(Wd/pointSampl));
JxSamp = zeros(floor(Ln/pointSampl),floor(Wd/pointSampl));
JySamp = zeros(floor(Ln/pointSampl),floor(Wd/pointSampl));

for ite=1:floor(Ln/pointSampl)
    for jd = 1:floor(Wd/pointSampl)
        eXSamp(ite,jd)=eXs(pointSampl*ite,pointSampl*jd);
        eYSamp(ite,jd)=eYs(pointSampl*ite,pointSampl*jd);
        JxSamp(ite,jd)=JxD(pointSampl*ite,pointSampl*jd);
        JySamp(ite,jd)=JyD(pointSampl*ite,pointSampl*jd);
    end
end

%% Plot Question 2 - Part(a)
% This section plots conductivity sigma(x,y), V(x,y), E(x,y), and J(x,y).

f3=figure;
figure(f3)
subplot(2,2,1)
surface(cMp2');
title('Conductivity');
subplot(2,2,2)
surface(V2Mp');
title('Voltage Vx,Vy');
subplot(2,2,3)
quiver(y,x,eXSamp,eYSamp,10);
title('Electric Field, Ex,Ey');
subplot(2,2,4)
quiver(y,x,JxSamp,JySamp,10);
title('Current Density, Jx,Jy');

%% Question 2 - Part (b)
% The function has been created to take inputs and meshing is also changed
% to check the effect on results. Mesh density is passed through the
% functions to plot results. 

% It is with the meshing size double of standard set.
Plotconductivity(40, 10e-2, 3);

% It is meshing size triple of standard set.
Plotconductivity(80, 10e-2, 3);

% It is meshing size half of standard set.
Plotconductivity(10, 10e-2, 3);

%% Question 2 - Part (c)
% Here, the bottle neck size and conductivity will be same but different dimmensions.
% The plot is done for different sizes of bottle neck such as normal, 
% longer and narrower bottle necks. This is done by passing the values to
% functions as done in the following. 

Plotconductivity(100, 10e-2, 3);

% It has longer bottle neck

Plotconductivity2(100, 10e-2, 3, 0.5, 0.375);

% It has narrower bottle neck
Plotconductivity2(100, 10e-2, 3, 0.3333, 0.5);

%% Question 2 - Part (d)
% Now, the results with differnce resitivity are compared inside of the
% box. The results obtained by applying difference method are  very close
% as required. The method has tested with different levels of restiveity
% such as higher and lower resistivity. The function is designed to show
% shorten code and re-use the code by calling function as done in the
% follwoing. 

% It is with higher resistivity as compared to the original
Plotconductivity(20, 10, 3);

% It is with lower resistivity as compared to the original
Plotconductivity(20, 10e-5, 3);