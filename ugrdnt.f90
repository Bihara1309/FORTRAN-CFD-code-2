!
! this subroutine calculates the velocity gradient
!
subroutine ugrdnt(ug12,twt,twb,utaut,utaub,u,un,us,&
                  yn,yp,dn,ds,ib,ie,id,den,mu,uold)
    real ug12(id),u(id),un(id),us(id)
    real yn(id),dn(id),ds(id),yp(id),uold(id)
    real ib,ie,twt,twb,utaut,utaub,mu
    real den,nu
    ibm1 = ib-1
    iem1 = ie-1    
    do 11 i=ibm1,iep1
        u(i) = uold(i)+0.5*(u(i)-uold(i))
 11 continue
 do 12 i=ibm1,iep1
     uold(i) = u(i)
     12 continue
    do 10 i=ibm1,iem1
        ug = (u(i+1)-u(i))/(dn(i)+ds(i+1))
        ug12(i) = ug
!        print*,ug12(i)
 10 continue
        ug12(ie) = (-u(ie)+u(ie+1))/(dn(ie)+ds(ie+1))
        twt = mu*abs(ug12(ie))
        utaut = sqrt(twt/den)
        twb = mu*abs(ug12(ib))
        utaub = sqrt(twb/den) 
!        print*,twb,twt
end subroutine ugrdnt