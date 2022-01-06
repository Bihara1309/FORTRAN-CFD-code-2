!
! this subroutine calculates the source term for the discrete equations
!
subroutine source(b,dpdx,yn,ys,dn,ds,ib,ie,id)
    real b(id),yn(id),ys(id)
    real ib,ie,dpdx
    real dn(id),ds(id)
    do 10 i=ib,ie
        vol = dn(i)+ds(i)
        b(i) = dpdx*vol
 10 continue
end subroutine source