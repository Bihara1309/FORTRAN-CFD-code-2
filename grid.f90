! this subroutine creates a non uniform grid with very fine mesh near the walls for a channel
!
subroutine grid(yp,yn,ys,ds,dn,ncv,h,ib,ie,id)
    real yp(id),yn(id),ys(id),ds(id)
    real dn(id)
    real ib,ie,ncv,h
    open (unit = 23, file = 'grid.dat', status = 'replace', form = 'formatted')
!
!
!    ib = 2
!    ncv = 40
!    h = 0.03
!    ie = ncv+1
    ibm1 = ib-1
    iep1 = ie+1
    ic = ib-1+((ie-ib+1)/2)
    icp1 = ic+1
    hh = h/2
    nn = ncv/2
    rlx = 1.2
    rk = hh/((rlx**nn)-1)
    yn(ibm1) = 0.0
    do 10 i=ib,ic
        yn(i) = ((rlx**(i-1))-1)*rk
        ys(i) = yn(i-1)
        dy = yn(i)-yn(i-1)
        dy2 = dy/2
        yp(i) = yn(i)-dy2
        ds(i) = dy2
        dn(i) = dy2
 10 continue
 ys(ibm1) = -yn(ib)
    do 11 i=icp1,ie
        j=ie-i+2
        yn(i) = yn(i-1)+(yn(j)-yn(j-1))
        ys(i) = yn(i-1)
        dy = yn(i)-yn(i-1)
        dy2 = dy/2
        yp(i) = yn(i)-dy2
        ds(i) = dy2
        dn(i) = dy2
 11 continue
 yn(ie) = h
 ys(iep1) = yn(ie)
        dn(ibm1) = dn(ib)
        yp(ibm1) = yn(ibm1)-dn(ibm1)
        ds(iep1) = dn(ie)
        yp(iep1) = yn(ie)+ds(iep1)
    write(23,200)
    do i=ibm1,iep1
        write(23,210) i,yp(i),yn(i),ys(i),dn(i),ds(i)
    end do
    200 format('',/,'grid parameters: ',//,5X,'    I  ',4X,'     yp   '&
               ,5X,'     yn   ',4X,'     ys    ',4X,'         dn    '&
               ,4X,'      ds    ')
    210 format('',5X,I5,5(3X,E14.7))

!    do 12 i=ibm1,iep1
!        print*,yp(i)
! 12 continue
end subroutine grid