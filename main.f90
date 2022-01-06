! this is the mainline code for the turbulent channel flow
program turbflow
    real ncv,id,ib,ie
    real mu,dpdx,den,h,kappa
    real an(100),as(100),ap(100),b(100)
    real yp(100),yn(100),ys(100),ug12(100)
    real u(100),un(100),us(100)
    real uold(100),ds(100),dn(100)
    real ypls(100),mutot(100),mut(100)
    real t(100),al(100),bt(100)
    real utaut,utaub,twt,twb,nu
    open (unit = 24, file = 'output.dat', status = 'replace', form = 'formatted')
!
!    
    ncv = 80
    h = 0.05
    ib = 2
    ie = ncv+1
    id = 100
    den = 1000
    kappa = 0.41
    mu = 0.001
    nu = mu/den
    dpdx = 15
    uinit = 0.1
    ibm1 = ib-1
    iep1 = ie+1
    do 11 i=ibm1,iep1
        u(i) = uinit
        uold(i) = uinit
 11 continue
!
! generate the non uniform grid
    call grid(yp,yn,ys,ds,dn,ncv,h,ib,ie,id)
!
!
! initialize the turbulence field
    call trbini(mutot,utaut,utaub,mut,ug12,ib,ie,id)
! start the solution loop
    do 10 i=1,20000
    call bcvel(ap,as,an,b,ib,ie,id)
    call trbvis(mu,mutot,mut,ypls,den,kappa,utaub,&
                utaut,ug12,yp,yn,ys,ds,dn,ncv,h,ib,ie,id)
!
    call coef(mutot,dn,ds,yp,yn,ys,an,as,ap,ib,ie,id,mu)
    call source(b,dpdx,yn,ys,dn,ds,ib,ie,id)
    call tdma(an,as,ap,b,u,al,bt,ib,ie,id)
    call ugrdnt(ug12,twt,twb,utaut,utaub,u,un,us,&
                yn,yp,dn,ds,ib,ie,id,den,mu,uold)
 10 continue
    write(24,100)
    do i=ib,ie
        write(24,110) i,yp(i),ypls(i),u(i),mutot(i),mut(i),ug12(i)
    end do
    100 format('',/,'Final Solution Fields: ',//,5X,'    I  ',4X,'      y  '&
               ,3X,'     y+n    ',5X,'      Velocity   ',1X,'  total viscosity    '&
               ,3X,' turb vis   ',4X,'     UG12   ')
    110 format('',5X,I5,6(3X,E14.7))
    write(24,120) utaub,utaut
    120 format('',/,'friction velocity (bottom wall):',E14.7,/,'',&
               'friction velocity (top wall):',E14.7)
end program turbflow