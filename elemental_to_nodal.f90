!Gets the nodal values from averaging values in each elemental
!That's a terrible explanation
!Currently only works for m=4, probably generalisable

subroutine elemental_to_nodal(elemental,nodal,inv_icon,np,nel,m)
integer np,nel,m
integer,dimension(m,np) :: inv_icon
real(8), dimension(np) :: nodal
real(8),dimension(nel) :: elemental
integer iel,ip,k

do ip=1,np

if (inv_icon(2,ip) .eq. 0) then !if node only adjacent to one element
nodal(ip) = elemental(inv_icon(1,ip)) !nodal value equal to value at that element
else if (inv_icon(3,ip) .eq. 0) then !if adjacent to two
nodal(ip) = (elemental(inv_icon(1,ip)) + elemental(inv_icon(2,ip)))/2.d0 !equal to average value of the two
else if (inv_icon(4,ip) .eq. 0) then !etc.
nodal(ip) = (elemental(inv_icon(1,ip)) + elemental(inv_icon(2,ip)) + elemental(inv_icon(3,ip)))/3.d0
else 
nodal(ip) = (elemental(inv_icon(1,ip)) + elemental(inv_icon(2,ip)) + elemental(inv_icon(3,ip)) + elemental(inv_icon(4,ip)))/4.d0
end if

end do

end subroutine