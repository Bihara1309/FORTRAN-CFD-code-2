!
! this subroutine initialize the turbulence field
subroutine trbini(mutot,utaut,utaub,mut,ug12,ib,ie,id)
    real mutot(id),mut(id),ug12(id)
    real utaut,utaub
    real ib,ie
    ibm1 = ib-1
    iep1 = ie+1
    utaut = 1.0
    utaub = 1.0
    do 10 i=ibm1,iep1
        mut(i) = 1.0
        mutot(i) = 2.0
        ug12(i) = 1.0
!        print*,mut(i)
 10 continue
end subroutine trbini