mass = 900.0;
Izz = 1100;
cog = [0.3462];
wheelbase = 2.6;
Cf = -62228.0;
Cr = -50997.0;
delta = deg2rad(10);
V = 35;

yawcompare (cog,wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr);

% timestep = 2;
% a = wheelbase*cog;
% b = wheelbase-a;
% [Ndelta,Ydelta,Nr,Yb,Nb,Yr] = derivatives(a,b,Cf,Cr);
% [k,wn,C2,Cc,c,z,wd,phi,X,tanp] = parameters(a,delta,V,Ndelta,Ydelta,Nr,Yb,Nb,Yr,mass,Izz,Cf,Cr);
% [r,time,rinf] = yawrate (a,b,timestep,delta,mass,Izz,V,Cf,Cr);
% plot (time,r)
% txt = string(rinf);
% text(1.5,rinf-0.05,txt)

%YAW RATE COMPARISON
function [] = yawcompare (cog,wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr)
labels = cell(1,length(cog))    
    for i = 1:1:length(cog)
        a = wheelbase*cog(i);
        b = wheelbase-a;
        [r,time,rinf,R,K,x,y,Vlim,vel,rd] = analysis(wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr);
        plot(time,r)
        txt = string(rinf);
        text(1.5,rinf-0.05,txt)
        labels{1,i} = string(cog(i))
        hold on
    end
    hold off
    legend(labels)
    legend('location','southeast')
end

function [] = velocitiescompare (cog,wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr)
labels = cell(1,length(cog))    
    for i = 1:1:length(cog)
        a = wheelbase*cog(i);
        b = wheelbase-a;
        [r,time,rinf,R,K,x,y,Vlim,vel,rd] = analysis(wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr);
        plot(vel,rd)
        title('Limit velocities')
        txt = string(Vlim);
        text(Vlim+0.05,rd(50),txt)
        xline(Vlim)
        labels{1,i} = string(cog(i))
        hold on
    end
    hold off
    legend(labels)
    legend('location','south')
end

function [] = radiuscompare (cog,wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr)
labels = cell(1,length(cog))
    for i = 1:1:length(cog)
        a = wheelbase*cog(i);
        b = wheelbase-a;
        [r,time,rinf,R,K,x,y,Vlim,vel,rd] = analysis(wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr);
        plot(x,y)
        title('Turning Trayectory')
        txt = string(R);
        text(R-1,R-1,txt)
        labels{1,i} = string(cog(i))
        hold on
    end
    hold off
    legend (labels)
    legend('location','south')
end

function[] = allcompare(cog,wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr)
    for i = 1:1:length(cog)
        a = wheelbase*cog(i);
        b = wheelbase-a;
        [r,time,rinf,R,K,x,y,Vlim,vel,rd] = analysis(wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr);
        figure(i);
        tiledlayout(3,1)
        nexttile
        plot(time,r)
        title('yaw rate response')
        txt = string(rinf);
        text(1.5,rinf-0.05,txt)

        nexttile
        plot(vel,rd)
        title('Limit velocities')
        txt = string(Vlim);
        text(Vlim+0.05,rd(50),txt)
        xline(Vlim)


        nexttile
        plot(x,y)
        title('Turning Trayectory')
        txt = string(R);
        text(1,1,txt)
        text(x(1),y(1),string(K))

end
end

%derivatives
function [Ndelta,Ydelta,Nr,Yb,Nb,Yr] = derivatives(a,b,Cf,Cr)
Ndelta = -a*Cf;
Ydelta = -Cf;
Nr = (a^2*Cf + b^2*Cr);
Yb = Cf+Cr;
Nb = a*Cf - b*Cr;
Yr = (a*Cf - b*Cr);
end

%parameters
function [k,wn,C2,Cc,c,z,wd,phi,X,tanp] = parameters(a,delta,V,Ndelta,Ydelta,Nr,Yb,Nb,Yr,mass,Izz,Cf,~)
k = Nb +(Yb*Nr/V)/(mass*V) - (Nb*Yr/V)/(mass*V);
wn = (k/Izz)^0.5;
C2 = (Ydelta*Nb)/(mass*V) - (Yb*Ndelta)/(mass*V);
Cc = 2*Izz*wn;
c = -(Nr/V + (Izz*Yb)/(mass*V));
z = c/Cc;
wd = ((1-z^2)^0.5)*wn;
tanp = -(-C2*wd*Izz)/(Cf*a*k + z*wn*C2*Izz);
phi = atan(tanp);
X = (-C2*delta)/(k*sin(phi));
end

%yawrate calculation
function [r,time,rinf] = yawrate (a,b,t,delta,mass,Izz,V,Cf,Cr)
[Ndelta,Ydelta,Nr,Yb,Nb,Yr] = derivatives(a,b,Cf,Cr);
[k,wn,C2,~,~,z,wd,phi,X,~] = parameters(a,delta,V,Ndelta,Ydelta,Nr,Yb,Nb,Yr,mass,Izz,Cf,Cr);
time = 0:0.02:t;
rinf = C2*delta/k;
if z < 1
    r = ((X*exp(-z*wn*time).*sin(wd*time+phi) + (C2/k)*delta));
end
if z == 1
    r = (((-C2*delta/k)+((-Cf*delta*a/Izz)-wn*(C2*d/k))*time)*exp(-wn*time)+C2*delta/k);
end
if z > 1
    f = (-z-(z^2-1)^0.5)*wn;
    g = (-z+(z^2-1)^0.5)*wn;
    A = 1/(g-f)*((-Cf*delta*a)/Izz + f*rinf);
    B = -(A+rinf);
    r = (A*exp(g*time)+B*exp(f*time)+rinf);
end
end

%corner radius
function [R,K,x,y] = cornerrad(mass,wheelbase,a,b,Cf,Cr,delta,V)
K = mass*9.81/wheelbase*(a/(Cr*pi/180) - b/(Cf*pi/180));
R = (1+K*V^2)/((delta/180*pi)/wheelbase);
angle = 0:0.01:90;
x = R*cos(deg2rad(angle));
y = R*sin(deg2rad(angle));
end

%speedlimits
function [Vlim,vel,rd] = speedlim(a,b,mass,wheelbase,Cf,Cr,K)
if (K > 0)
    Vlim = (1/(K/180*pi))^0.5;
end
if(K < 0)
    Vlim = (-1/(K/180*pi))^0.5;
end
vel = 0:1:100;
rd = (vel/wheelbase)./(1+(K/180*pi)*vel.^2);
end
    
        
function [r,time,rinf,R,K,x,y,Vlim,vel,rd] = analysis(wheelbase,a,b,timestep,delta,mass,Izz,V,Cf,Cr,frontw)
[r,time,rinf] = yawrate (a,b,timestep,delta,mass,Izz,V,Cf,Cr);
[R,K,x,y] = cornerrad(mass,wheelbase,a,b,Cf,Cr,V,delta);
[Vlim,vel,rd] = speedlim(a,b,mass,wheelbase,Cf,Cr,K);

end



