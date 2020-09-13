clear all ;clc;close all;

%% Engine specification
b = 0.137 ;% Bore in meter
l = 0.165 ;%Stroke Length in meter
c_r = 14.5;%Compression ratio
d_0 = 0.15*10^(-3);%Injector hole diameter in meter
alpha = 50;% Angle between spray centreline and cylinder head in degrees
beta  = 45;%Half slice angle
b_d = 0.093;%Bowl diameter in meter
b_h =0.017;%Depht of the bowl in meter
rho_a = 830;%Density of incylinder gas at start of ignition
delta_p = 125*10^5;%Pressure difference between injector and ambient
mu_a = 2.98*10^(-6)*rho_a;%Dynamics viscosity of incylinder gas at start of ignition

%% Prelimanry Calculations
v_s = pi*((b/2)^2)*l;%Swept volume
v_c = v_s/(c_r-1);%Clearance volume
Vcyl = 1.3*v_c; % Instantaneous volume of cylinder
hp_2 = Vcyl/(pi*(b/2)^2);%Distance between piston top and cylinder head
hp_0 = pi*((b_d/2)^2)*b_h;%Distance betweem the piston top and piston bowl
x_t = (hp_2*tand(alpha)+(b/2))*cosd(alpha);%Maximum distance upto which plane C can move
hp_1 = hp_2 + hp_0 ;
dx_t = x_t/3;

x = 0.01;
%for x = 0:dx_t:x_t
    %% Finding the coordinates of the A,B,C,D,E,O
    A = [(x/sind(alpha)-hp_1)*tan(alpha),-hp_1];
    B = [(x/sind(alpha)-hp_2)*tan(alpha),-hp_2];
    C = [x/cos(alpha),0];
    D = [b/2,((b/2)/tand(alpha))-(x/sind(alpha))];
    E = [0,-x/sind(alpha)];
    O = [x*cosd(alpha),-x*sind(alpha)];
    
    %% Finding lengths of OA,OB,OC,OD,OE
    OA = distance(O,A);
    OB = distance(O,B);
    OC = distance(O,C);
    OD = distance(O,D);
    OE = distance(O,E);
    
    %% Finding l1,l2,h1,h2
    l_1 = 0;
    l_2 = 0;
    rp_0 = b_d/2;
    switch l_1
        case D(2) > C(2)
            l_1 = abs(OC);
        case (O(2) <= D(2) & D(2) < C(2))
            l_1 = abs(OD);
        case D(2) < O(2)
            l_1 = -abs(OD);
    end
    
    switch l_2
        case A(1) <= 0
            l_2 = abs(OE);
        case (B(1) >= rp_0 & O(2)>= B(2))
            l_2 = abs(OB);
        case B(1) >= rp_0 & O(2) <B(2)
            l_2 = -abs(OB);
        case ( A(1) > 0 & B(1) < rp_0 & O(2) >= A(2) )
            l_2 = abs(OA);
        case ( A(1) > 0 & B(1) < rp_0 & O(2)< A(2))
            l_2 = -abs(OA);
    end
    
    s_1 = x/(sind(alpha)*cosd(alpha));
    s_2 = x*tand(alpha);
    alpha_c1 = atan2(sind(alpha)*tand(alpha),1); % Angle which determinaes the shape of traingular boundary
    h_1 = (l_1+s_1-s_2)*tan(alpha_c1);
    h_2 = (s_1-s_2-l_2)*tan(alpha_c1);
    
    theta = 0.025*((rho_a*delta_p*d_0^2)/(mu_a^2))^0.25;
    t = 3;
    rho_0 = 1.225;%Density of air at SATP in kg/m^2
    Lp_0 = 4.347*((rho_a/rho_0)^(-0.175))*((2*delta_p/rho_a)^0.25)*sqrt(d_0*t);%Spray tip peneration at SOC
    x_0 = Lp_0/(1+tan(theta)); %Height of spray cone
    r_0 = x_0*tan(theta);%Diameter of hemisphere
    phi = 0.7;%Equivalence ratio
    SL_0 = 0.4;%laminar flame speed of NG-Air mixture
    
    %% Calculation of r_x
    Lp = Lp_0 + SL_0*t;%Distance between injector hole and furthest location of flame front raltvie to injector
    r = Lp - x_0;%Radius of hemisphere front at any instant t
    r_x = 0;%Distance between piston top and cylinder head
    switch r_x
        case (0 < x) & (x <= x_0)
            r_x = x*tan(theta)+(r-r_0);
        case (x_0 <x) & (x <= (x_0+r))
            r_x = sqrt(r^2 + (x - x_0)^2);
        case x > (x_0 +r)
            r_x = 0;
    end
    
    %% Calculaiton of l_x
    dx_c = r_x/10;
    l_x = 0;
    for x_c = -r_x:dx_c:r_x
       %% Finding y_l and y_c  
       y_l = 0;
       y_c = 0;
       switch y_l
           case x_c <= -l_1
               y_l = 0;
           case -l_1 < x_c & x_c <= l_2
               y_l = (x_c - l_2)*(h_2-h_1)/(l_1+l_2) + h_2;
           case x_c > l_2
               y_l = 0;
       end
       
       switch y_c
           case x_c <= -r_x
               y_c = 0;
           case (-r_x < x_c) & (x_c <= r_x)
               y_c = sqrt(r_x^2 - x_c^2);
           case x_c > r_x
               y_c = 0;
       end
       
       syms p
       if (y_l >= y_c)
           f = (2*r_x)/sqrt(r_x^2 - p^2);
       else
           f = 0;
       end
       l_x =  simplify(int(f,p,-r_x,r_x));
    end
    
    %% Finding the total flame area
    Area = 0;
    switch Area
        case  x <= x_0
            Area = (2*l_x*r*x_t)/cos(theta);
        case  x > x_0
            Area = (2*l_x*r*x_t)/r_x;
    end
    
   
%end
j = 0:0.001:pi;% Theta for the semicircle coordinates
plot(x_c*cos(j),x_c*sin(j))%Semicircle
hold on
plot([l_2,l_2],[0,h_2])%Right vertical constrain
plot([-l_1,-l_1],[0,h_1])%Left vertical constrain
plot([l_2,-l_1],[h_2,h_1])%Inclined constrain







