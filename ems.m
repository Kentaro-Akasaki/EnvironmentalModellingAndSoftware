clear

%Import the observed data
tbl_layer1=readtable("data_layer1.csv");
tbl=readtable("data_20170822.csv");
syms x y z k2(x) k2t(x) p1(y) f(x) j0(x) n0(x) p2(x) p3(z) const(s);

% Estimation of the classical diffusion coefficient of layer 1
x0=824.1;
T=0.0921;
lambda=0.001014;
x1=668.8;
l1=5.5;
myfittype=fittype(@(a,x) T*(x0*cos(a*x)+((x1-x0*cos(a*l1))/sin(a*l1))*sin(a*x)));

fit_result=fit(tbl_layer1.depth,tbl_layer1.concentration,myfittype,"StartPoint",[0.1]);
a=coeffvalues(fit_result);
k1=lambda/a^2;

% Setting the values of the classical diffusion coefficients
k3=1.386;
at=0.22;
bt=-1.144;
k2(x)=at*x+bt;

%Setting the values of each parameters
x1=668.8;
x2=449.5;
T=0.0921; %The value of T part at t=2352[day]

l1=5.5;
l2=11.5;
L=17.5;

%Analytic solution of the layer 1
p1(y)=T*(x0*cos(sqrt(lambda/k1)*y)+((x1-x0*cos(sqrt(lambda/k1)*l1))/sin(sqrt(lambda/k1)*l1))*sin(sqrt(lambda/k1)*y));

%Analytic solution of the layer 2
f(x)=2*sqrt(lambda)*sqrt(k2)/at;
j0(x)=besselj(0,f(x));
n0(x)=bessely(0,f(x));
c2=(x2*j0(l1)-x1*j0(l2))/(n0(l2)*j0(l1)-j0(l2)*n0(l1));
c1=(x1-c2*n0(l1))/j0(l1);

p2(x)=T*(c1*j0(x)+c2*n0(x));

%Analytic solution of the layer 3
h=sqrt(lambda/k3);
u=x2/(cos(h*l2)+tan(h*L)*sin(h*l2));
v=u*tan(h*L);
p3(z)=T*(u*cos(h*z)+v*sin(h*z));

%The solution with constant coefficient [Akasaki, et al., 2022]
const(s)=0.0226*2622*cos(0.045*s);

%Creating the figure of the observed data and the analytic solutions
fig=figure('Position', [100 100 1280 720]);
sc=scatter(tbl.depth,tbl.concentration,"filled","LineWidth",3);
sc.SizeData=100;
hold on
fplot(p1,[0 5.5],"LineWidth",4)
fplot(p2,[5.5 11.5],"LineWidth",4)
fplot(p3,[11.5 17.5],"LineWidth",4)
fplot(const,[0 17.5],"--","LineWidth",4)
legend1=legend("Observed data","p_1","p_2","p_3","Constant coef.");
set(legend1,'FontSize',25,'FontWeight','bold');
xlabel("Distance from bottom sediments [m]",'FontSize',27,'FontWeight','bold');
ylabel("Activity concentration of ^{137}Cs [Bq/m^3]",'FontSize',27,'FontWeight','bold');

box on;
ax=gca;
ax.LineWidth=3;
ax.FontSize=25;
ax.FontWeight="bold";
hold off
fig.PaperOrientation='Landscape';
