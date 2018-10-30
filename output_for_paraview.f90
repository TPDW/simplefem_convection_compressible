subroutine output_for_paraview (np,nel,x,y,Vx,Vy,p,T,density,qx,qy,icon,istep,dxu_nodal,dyv_nodal,dxv_nodal,dyu_nodal,dxu_elemental,dyv_elemental,dxv_elemental,dyu_elemental,phi,phi_el,adiabatic)
implicit none
integer np,nel
real(8), dimension(np)    :: x,y,Vx,Vy,T,density
real(8), dimension(np)    :: dxu_nodal,dyv_nodal,dxv_nodal,dyu_nodal,phi
real(8), dimension(nel)   :: p,qx,qy
real(8), dimension(nel)   :: dxu_elemental,dyv_elemental,dxv_elemental,dyu_elemental,phi_el,adiabatic
integer, dimension(4,nel) :: icon
integer istep
character(len=6) cistep

integer i,iel

!=======================================

call int_to_char(cistep,6,istep)

open(unit=123,file='OUT/Paraview/solution_'//cistep//'.vtu',status='replace',form='formatted')
write(123,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(123,*) '<UnstructuredGrid>'
write(123,*) '<Piece NumberOfPoints="',np,'" NumberOfCells="',nel,'">'
!.............................
write(123,*) '<PointData Scalars="scalars">'

write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="Velocity" Format="ascii">'
do i=1,np
write(123,*) Vx(i),Vy(i),0
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" Name="density" Format="ascii">'
do i=1,np
write(123,*) density(i)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" Name="Temperature" Format="ascii">'
do i=1,np
write(123,*) T(i)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" NumberOfComponents="4" Name="Strain" Format="ascii">'
do i=1,np
write(123,*) dxu_nodal(i),dyv_nodal(i),dxv_nodal(i),dyu_nodal(i)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" Name="phi" Format="ascii">'
do i=1,np
write(123,*) phi(i)
end do
write(123,*) '</DataArray>'

write(123,*) '</PointData>'
write(123,*) '<Points>'
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
do i=1,np
write(123,*) x(i),y(i),0.d0
end do
write(123,*) '</DataArray>'
write(123,*) '</Points>'
!.............................
write(123,*) '<CellData Scalars="scalars">'

write(123,*) '<DataArray type="Float32" Name="pressure" Format="ascii">'
do iel=1,nel
write(123,*) p(iel)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" Name="adiabatic heating" Format="ascii">'
do iel=1,nel
write(123,*) adiabatic(iel)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="heat flux" Format="ascii">'
do iel=1,nel
write(123,*) qx(iel),qy(iel),0
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" NumberOfComponents="4" Name="Elemental Strain" Format="ascii">'
do i=1,nel
write(123,*) dxu_elemental(i),dyv_elemental(i),dxv_elemental(i),dyu_elemental(i)
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Float32" Name="phi elemental" Format="ascii">'
do iel=1,nel
write(123,*) phi_el(iel)
end do
write(123,*) '</DataArray>'


write(123,*) '</CellData>'
!.............................
write(123,*) '<Cells>'
write(123,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do iel=1,nel
write(123,*) icon(1,iel)-1,icon(2,iel)-1,icon(3,iel)-1,icon(4,iel)-1
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
write(123,*) (iel*4,iel=1,nel)
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Int32" Name="types" Format="ascii">'
write(123,*) (9,iel=1,nel)
write(123,*) '</DataArray>'
write(123,*) '</Cells>'
write(123,*) '</Piece>'
write(123,*) '</UnstructuredGrid>'
write(123,*) '</VTKFile>'
close(123)


end subroutine
