! this subroutine evaluates the coefficients of the discrete equations
!
subroutine coef(mutot,dn,ds,yp,yn,ys,an,as,ap,ib,ie,id,mu)
    real mutot(id)
    real an(id),as(id),ap(id)
    real dn(id),ds(id),yp(id),yn(id),ys(id)
    real ib,ie,mu
!    an(ib) = mutot(ib)/(yn(ib)-yp(ib))
!    as(ib) = mu/(yp(ib)-ys(ib))
!    ap(ib) = an(ib)+as(ib)
    do 10 i=ib,ie
        an(i)=mutot(i)/(yp(i+1)-yp(i))
        as(i)=mutot(i-1)/(yp(i)-yp(i-1))
        ap(i)=an(i)+as(i)
!        print*,i,an(i),as(i),ap(i)
 10 continue
! an(ie) = mu/dn(ie)
! as(ie) = mutot(ie-1)/ds(i)
! ap(ie) = an(ie)+as(ie)
end subroutine coef