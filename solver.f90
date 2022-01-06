!
! this subroutine is a Gauss Siedel solver for the turbulent channel flow\
!
subroutine solver(an,as,ap,b,u,un,us,uold,yn,ib,ie,id)
    real u(id),un(id),us(id),uold(id),b(id)
    real an(id),as(id),ap(id),yn(id)
    real ib,ie
    ibm1 = ib-1
    iep1 = ie+1
    u(ibm1) = an(ibm1)*uold(ib)
!    print*,u(ibm1)
    do 10 i=ib,ie
        u(i)=uold(i)+&
        (0.5*(((1/ap(i))*((an(i)*uold(i+1))+(as(i)*u(i-1))+b(i)))-uold(i)))
!        print*,u(i)
 10 continue
        u(iep1) = as(iep1)*u(ie)
!        print*,u(iep1)
    do 11 i=ibm1,iep1
        uold(i)=u(i)
 11 continue
! print*,u(ie)
end subroutine solver