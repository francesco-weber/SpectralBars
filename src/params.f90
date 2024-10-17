    subroutine param
    use mod_main
    implicit none
    integer :: i, j
    real*8 :: rs
  
    pi = acos(-1.d0)

    
    zi=(zero,1.d0)
    do i=0,nxm
        zilx(i) = zi*i*rlx
    enddo
    do j=0,nym
        zily(j) = 0.5d0*zi*j*pi
    enddo
    dslog = 2.5d0*log(2.5d0*ds)
    c0 = 1.d0/(6.d0-dslog)**(2.d0)
    fs2 = teta0/c0
    rs = 1.65d0
    tetacr = 0.047d0
    f02 = rs*fs2*ds
    q0 = ds/(0.4d0*sqrt(fs2))
    rsb = rt/beta
   
    return
    end