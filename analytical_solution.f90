
!pth = -(1.d0/x**2+1.d0/y**2)-x*y



!Cases
!1 - 2D Cartesian Linear
!2 - 2D Cartesian Sinusoidal
!3 - 1D Cartesian Linear
!4 - Arie van den Berg
!5 - 1D Cartesian Sinusoidal


function uth (x,y,ibench)
implicit none
real(8) uth,x,y,r,phi
integer ibench
r=sqrt(x**2+y**2)
phi = atan(y/x)
select case(ibench)
case(1)
uth = 1.d0/x 
case(2)
uth = 1.d0/cos(x)
case(3)
uth = 1.d0/x
case(4)
uth=1.d0-x
case(5)
uth = 1.d0/cos(x)
end select
end function

function vth (x,y,ibench)
implicit none
real(8) vth,x,y,r,phi
integer ibench
select case(ibench)
case(1)
vth = 1.d0/y
case(2)
vth = 1.d0/cos(y)
case(3)
vth=0.d0
case(4)
vth=0.d0
case(5)
vth=0.d0

end select
end function

function pth (x,y,ibench)
implicit none
real(8) pth,x,y,r,phi
integer ibench
select case(ibench)
case(1)
pth = -(1.d0/x**2+1.d0/y**2)-x*y + 3.25
case(2)
pth = (sin(x)/cos(x)**2+sin(y)/cos(y)**2) - (sin(x) + sin(y)) - 0.78223604709813
case(3)
pth =  -1.d0/x**2 - x**2/2.d0 + 5.d0/3.d0
case(4)
pth = 10.d0*log(x-1) -29.9 
!pth = 1.d0/(1.d0-x)*x
case(5)
pth = sin(x)/cos(x)**2-sin(x) - 0.391118
end select

end function

function rho(x,y,ibench)
implicit none
real(8) rho,x,y,r,phi
integer ibench
select case(ibench)
case(1)
rho=x*y
case(2)
rho = cos(x)*cos(y)
case(3)
rho=x
case(4)
rho = 1.d0/(1.d0-x)
case(5)
rho=cos(x)
end select
end function

function drhodx(x,y,ibench)
implicit none
real(8) drhodx,x,y,r,phi
integer ibench
select case(ibench)
case(1)
drhodx=y
case(2)
drhodx=-sin(x)*cos(y)
case(3)
drhodx=1.d0
case(4)
drhodx=1.d0/(1.d0-x)**2
case(5)
drhodx=-sin(x)
end select
end function


function drhody(x,y,ibench)
implicit none
real(8) drhody,x,y,r,phi
integer ibench
select case(ibench)
case(1)
drhody=x
case(2)
drhody=-cos(x)*sin(y)
case(3)
drhody=0.d0
case(4)
drhody=0.d0
case(5)
drhody=0.d0
end select
end function

function gx(x,y,ibench)
implicit none
real(8) gx,x,y
integer ibench
select case(ibench)
case(1)
gx=-1.d0/x
case(2)
gx = -1.d0/cos(y)
case(3)
gx=-1.d0
case(4)
gx=10.d0
case(5)
gx=-1.d0
end select
end function

function gy(x,y,ibench)
implicit none
real(8) gy,x,y
integer ibench
select case(ibench)
case(1)
gy=-1.d0/y
case(2)
gy=-1.d0/cos(x)
case(3)
gy=0.d0
case(4)
gy=0.d0
case(5)
gy=0.d0
end select
end function