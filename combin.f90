module combin
implicit none
contains
subroutine ranper(a,setup)
! Random permutation code 
! adapted from Nijenhuis & Wilf
integer, dimension(:), intent(inout):: a
logical, intent(in):: setup
real,allocatable::xrand(:)
integer :: m,n,l,l1
n=size(a)
if (setup) then
  do m=1,n
    a(m)=m
  enddo
endif
allocate (xrand(n))
call random_number(xrand(:))
do m=1,n
  l=m+floor(xrand(m)*(n+1-m))
  l1=a(l)
  a(l)=a(m)
  a(m)=l1
enddo
deallocate (xrand) 
return
end subroutine ranper
end module combin
