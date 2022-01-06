!
! this subroutine calculates the turbulent viscosity of the channel flow
subroutine trbvis(mu,mutot,mut,ypls,den,kappa,utaub,&
                  utaut,ug12,yp,yn,ys,ds,dn,ncv,h,ib,ie,id)
!
    real mu,kappa,den,ncv,h
    real dn(id),ds(id),ug12(id)
    real yp(id),yn(id),ys(id)
    real mutot(id),ypls(id),mut(id)
    real utaub,utaut
    real ib,ie
    a = 25.0
    nu = mu/den
    ic = ib-1+((ie-ib+1)/2)
    icp1 = ic+1
    ibm1 = ib-1
    iep1 = ie+1
!
    do 10 i=ib,ic
        ypls(i) = (den*yn(i)*utaub)/mu
        damp = (1-exp(-ypls(i)/a))*kappa*yn(i)
        damp2 = damp**2
        mut(i) = den*damp2*abs(ug12(i))
        mutot(i) = mu+mut(i)
!        print*, mutot(i)
!        print*,ypls(i),mutot(i),dn(i)
 10 continue
! mutot(ic) = mu
 mutot(ibm1) = mutot(ib)
    do 11 i=icp1,ie
        ypls(i) = (den*(h-yn(i))*utaut)/mu
        damp = (1-exp(-ypls(i)/a))*kappa*(h-yn(i))
        damp2 = damp**2
        mut(i) = den*damp2*abs(ug12(i))
        mutot(i) = mu+mut(i)
!        print*, mutot(i)
 11 continue
 mutot(iep1) = mutot(ie)
! print*,utaub,utaut
end subroutine trbvis