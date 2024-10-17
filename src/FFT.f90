subroutine prefft(wd,wo,n)
    use mod_main
    implicit none
    
    integer,intent(in) :: n
    complex*16,intent(out) :: wd(n/2-1),wo(n/2-1)
    integer*2 :: irv
    integer :: i, ig, ipo, ipd
    real*8 :: argo, argd
      
    ig=nint(log(float(n))/log2_v)
  
    do i=1,n/2-1
        ipo = irv(i,ig-1)
        ipd = i
        argo = 2.d0*ipo*pi/n
        argd = 2.d0*ipd*pi/n
        wo(i) = cmplx(cos(argo),sin(argo))
        wd(i) = cmplx(cos(argd),-sin(argd))
    enddo
  
    return
    end
  
    ! **********************************************************************
  
    subroutine prefx(wx,wy,n)
    use mod_main
    implicit none

    integer,intent(in) :: n
    complex*16,intent(out) :: wx(n/2-1),wy(n/2-1)
    integer :: i
    real*8 :: arg

    do i=1,n/2-1
        arg = i*pi/n
        wx(i) = cmplx(sin(arg),cos(arg))
        wy(i) = cmplx(-sin(arg),cos(arg))
    enddo
  
    return
    end
    
    ! **********************************************************************
    
    integer*2 function irv(iw,nbit)
    implicit none
    integer,intent(in) :: iw,nbit
    integer :: i, iw1, iw2
  
    iw1 = iw
    irv = 0
    do i=1,nbit
        iw2 = iw1/2
        irv = 2*irv+(iw1-iw2*2)
        iw1 = iw2
    enddo
  
    return
    end

    ! *********************************************************************
    ! (used in atruod and atrvod for post-processing)
    subroutine ord(x,n,ig)
    use mod_main, only: ijd
    implicit none
  
    integer,intent(in) :: n, ig
    complex*16,intent(out) :: x(0:ijd-1)
    complex*16 :: zt
    integer :: i, k
    integer*2,external :: irv
  
    do i=1,n-2
        k = irv(i,ig)
        if(k.gt.i)then
            zt = x(k)
            x(k) = x(i)
            x(i) = zt
        endif
    enddo
  
    end
    
    ! *********************************************************************

    subroutine xtof(f,x,n,w)
    use mod_main, only: zero
    implicit none
  
    integer,intent(in) :: n
    complex*16,intent(out) :: f(0:n-1)          !mt CHECK DIMENSIONS!!!
    complex*16,intent(in) :: x(0:n/2-1), w(n/2-1)    !mt CHECK DIMENSIONS!!!
    complex*16 :: zt
    integer :: k
  
    f(0) = cmplx(real(x(0)),real(x(0)))
    f(n/2) = (zero,zero)
    !!$omp parallel do private(k) 
    do k=1,n/2-1
        zt = x(k)*w(k)
        f(k) = x(k)+zt                    ! first half (from 1 to n/2-1)
        f(n-k) = conjg(x(k))-conjg(zt)    ! second half (from n/2+1 to n-1)
    enddo
  !!$omp end parallel do
    return
    end
  
    ! *********************************************************************
  
    subroutine ftox(x,f,n,w)
    implicit none
 
    integer,intent(in) :: n
    complex*16,intent(in) :: f(0:n-1), w(n/2-1)    !mt CHECK DIMENSIONS!!!
    complex*16,intent(out) :: x(0:n/2-1)   !mt CHECK DIMENSIONS!!!
    integer :: k
  
    x(0) = real(f(0))
    !!$omp parallel do private(k) 
    do k=1,n/2-1
        x(k) = 0.25d0*(f(k)+conjg(f(n-k))-(f(k)-conjg(f(n-k)))*w(k))
    enddo
    !!$omp end parallel do
    
    return
    end
    
    ! **********************************************************************
  
    subroutine fftaod(x,n,ig,w)
    implicit none
    
    integer,intent(in) :: n, ig
    complex*16,intent(in) :: w(n/2-1)   !mt CHECK DIMENSIONS!!!
    complex*16,intent(inout) :: x(0:n-1)   !mt CHECK DIMENSIONS!!!
    integer :: n2, kl, kh, l, nb, ib, is, k
    complex*16 :: xx, wib
  
    n2=n/2
    ! l=1     nb=1
    do kl=0,n2-1
        kh=kl+n2
        xx=x(kh)
        x(kh)=x(kl)-xx
        x(kl)=x(kl)+xx
    enddo
    ! l=2,ig
    nb=2
    n2=n2/2
    do l=2,ig
        ! ib=0
        do kl=0,n2-1
            kh=kl+n2
            xx=x(kh)
            x(kh)=x(kl)-xx
            x(kl)=x(kl)+xx
        enddo
        ! ib=1,nb-1
        do ib=1,nb-1
            is=2*ib*n2
            wib=w(ib)
            do k=0,n2-1
                kl=k+is
                kh=kl+n2
                xx=wib*x(kh)
                x(kh)=x(kl)-xx
                x(kl)=x(kl)+xx
            enddo
        enddo
        nb=nb*2
        n2=n2/2
    enddo
  
    return
    end

    ! **********************************************************************
  
    subroutine ffttdo(x,n,ig,w)
    implicit none
  
    integer,intent(in) :: n, ig
    complex*16,intent(in) :: w(n/2-1)   !mt CHECK DIMENSIONS!!!
    complex*16,intent(inout) :: x(0:n-1)   !mt CHECK DIMENSIONS!!!
    
    integer :: n2, kl, kh, l, nb, ib, is, k
    complex*16 :: xx, wib
    real*8 :: fl
    
    n2=n/2
    !     l=1     k=0     nb=1
    do ib=0,n2-1
        kl=2*ib
        kh=kl+1
        xx=x(kh)
        x(kh)=x(kl)-xx
        x(kl)=x(kl)+xx
    enddo
    !     l=2,ig
    nb=2
    n2=n2/2
    do l=2,ig
        do ib=0,n2-1
            is=2*ib*nb
    !     k=0
            kl=is
            kh=kl+nb
            xx=x(kh)
            x(kh)=x(kl)-xx
            x(kl)=x(kl)+xx
    !     k=1,nb-1
            do k=1,nb-1
                kl=k+is
                kh=kl+nb
                xx=w(k*n2)*x(kh)
                x(kh)=x(kl)-xx
                x(kl)=x(kl)+xx
            enddo
        enddo
        nb=nb*2
        n2=n2/2
    enddo
  
    fl=1./n
    do k=0,n-1
        x(k)=x(k)*fl
    enddo
 
    return
    end

    
! **********************************************************************
! subroutines for computing transform (tr...) and "anti"-transform (atr...)
! **********************************************************************

    ! **********************************************************************
    subroutine atruod(z,zr,wxfx,wxfy,wox,woy,nx,ny)
    use mod_main, only: id, jd, ijd, idm, jdm, ijdm, zero, post_map,igx, igy, nxm, nym
    implicit none
  
    complex*16 :: z(0:idm,0:jdm)
    complex*16 :: wxfx(idm), wxfy(jdm), wox(idm), woy(jdm)
    complex*16 :: f(0:ijd-1)
    complex*16 :: x(0:ijdm)
    real*8 :: zr(0:2*id-1,0:2*jd-1)
    
    integer :: nx, ny
    integer :: i, j
    !integer :: i, j, igx, igy, nxm, nym
  
    !igx = nint(log(float(nx))/log(2.d0))
    !igy = nint(log(float(ny))/log(2.d0))
    !nxm = nx/2-1
    !nym = ny/2-1
  
    !write(*,*) nx !mt
    !$omp parallel do private(i, x, f) shared(z, zr, wxfx, wxfy, wox, woy, nx, ny) 
    !!$omp parallel do private(x) shared(zr, z) 
    do j=0,nym
        if(mod(j,2).eq.0) then
            do i=0,nxm
                x(i) = z(i,j)
            enddo
        else
            do i=0,nxm
                x(i) = cmplx(imag(z(i,j)),-real(z(i,j)))
            enddo
        endif
        call xtof (f,x,nx,wxfx)
        call fftaod(f,nx,igx,wox)
        if(post_map) then
            call ord(f,nx,igx)
        endif
        do i=0,nx-1
            zr(2*i,j) = real(f(i))
            zr(2*i+1,j) = imag(f(i))
        enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(j, x, f) shared(z, zr, wxfx, wxfy, wox, woy, nx, ny)
    do i=0,2*nx-1
        do j=0,nym
            if(mod(j,2).eq.0) then
                x(j) = zr(i,j)
            else
                x(j) = cmplx(zero,zr(i,j))
            endif
        enddo
        call xtof(f,x,ny,wxfy)
        call fftaod(f,ny,igy,woy)
        if(post_map) then
            call ord(f,ny,igy)
        endif
        do j=0,ny-1
            zr(i,2*j) = real(f(j))
            zr(i,2*j+1) = imag(f(j))
        enddo
    enddo
    !$omp end parallel do

    return
    end
  
    ! **********************************************************************
   
    subroutine traudo(zr,z,wfxx,wfxy,wdx,wdy,nx,ny)
    use mod_main, only: id, jd, ijd, idm, jdm, ijdm, igx, igy, nxm, nym
    implicit none
  
    integer :: nx, ny
    integer :: i, j

    complex*16 z(0:idm,0:jdm)
    complex*16 wfxx(idm),wfxy(jdm),wdx(idm),wdy(jdm)
    complex*16 f(0:ijd-1)
    complex*16 x(0:ijdm)
    real*8 zr(0:2*id-1,0:2*jd-1)
  
    !igx = nint(log(float(nx))/log(2.))
    !igy = nint(log(float(ny))/log(2.))
    !nxm = nx/2-1
    !nym = ny/2-1
  !$omp parallel do private(j, f, x) shared(zr, wfxx, wfxy, wdx, wdy, nx, ny)
    do i=0,2*nx-1
        do j=0,ny-1
            f(j) = cmplx(zr(i,2*j),zr(i,2*j+1))
        enddo
        call ffttdo(f,ny,igy,wdy)
        call ftox(x,f,ny,wfxy)
        do j=0,nym
            if(mod(j,2).eq.0)then
                zr(i,j) = real(x(j))
            else
                zr(i,j) = imag(x(j))
            endif
        enddo
    enddo
  !$omp end parallel do

  !$omp parallel do private(i, f, x) shared(zr, wfxx, wfxy, wdx, wdy, nx, ny)
    do j=0,nym
        do i=0,nx-1
            f(i) = cmplx(zr(2*i,j),zr(2*i+1,j))
        enddo
        call ffttdo(f,nx,igx,wdx)
        call ftox (x,f,nx,wfxx)
        if(mod(j,2).eq.0)then
            do i=0,nxm
                z(i,j) = x(i)
            enddo
        else
            do i=0,nxm
                z(i,j) = cmplx(-imag(x(i)),real(x(i)))
            enddo
        endif
    enddo
  !$omp end parallel do
    return
    end
  
    ! **********************************************************************
 
    subroutine atrvod(z,zr,wxfx,wxfy,wox,woy,nx,ny)
    use mod_main, only: id, jd, ijd, idm, jdm, ijdm, zero, post_map, igx, igy, nxm, nym
    implicit none
  
    integer :: nx, ny
    integer :: i, j
  
    complex*16 z(0:idm,0:jdm)
    complex*16 wxfx(idm),wxfy(jdm),wox(idm),woy(jdm)
    complex*16 f(0:ijd-1)
    complex*16 x(0:ijdm)
    real*8 zr(0:2*id-1,0:2*jd-1)
  
    !igx = nint(log(float(nx))/log(2.d0))
    !igy = nint(log(float(ny))/log(2.d0))
    !nxm = nx/2-1
    !nym = ny/2-1

    !$omp parallel do private(i, x, f) shared(z, zr, wxfx, wxfy, wox, woy, nx, ny)
    do j=0,nym
        if(mod(j,2).ne.0)then
            do i=0,nxm
                x(i) = z(i,j)
            enddo
        else
            do i=0,nxm
                x(i) = cmplx(imag(z(i,j)),-real(z(i,j)))
            enddo
        endif
        call xtof (f,x,nx,wxfx)
        call fftaod(f,nx,igx,wox)
        if(post_map) then
            call ord(f,nx,igx)
        endif
        do i=0,nx-1
            zr(2*i,j) = real(f(i))
            zr(2*i+1,j) = imag(f(i))
        enddo
    enddo
    !$omp end parallel do
  
    !$omp parallel do private(j, x, f) shared(z, zr, wxfx, wxfy, wox, woy, nx, ny)
    do i=0,2*nx-1
        do j=0,nym
            if(mod(j,2).ne.0) then
                x(j) = zr(i,j)
            else
                x(j) = cmplx(zero,zr(i,j))
            endif
        enddo
        call xtof(f,x,ny,wxfy)
        call fftaod(f,ny,igy,woy)
        if(post_map) then
            call ord(f,ny,igy)
        endif
        do j=0,ny-1
            zr(i,2*j) = real(f(j))
            zr(i,2*j+1) = imag(f(j))
        enddo
    enddo
  !$omp end parallel do

    return
    end
 
    ! **********************************************************************
  
    subroutine travdo(zr,z,wfxx,wfxy,wdx,wdy,nx,ny)
    use mod_main, only: id, jd, ijd, idm, jdm, ijdm, igx, igy, nxm, nym
    implicit none
  
    integer :: nx, ny
    !integer :: i, j, igx, igy, nxm, nym
    integer :: i, j
  
    complex*16 z(0:idm,0:jdm)
    complex*16 wfxx(idm),wfxy(jdm),wdx(idm),wdy(jdm)
    complex*16 f(0:ijd-1)
    complex*16 x(0:ijdm)
    real*8 zr(0:2*id-1,0:2*jd-1)
  
    !igx = nint(log(float(nx))/log(2.))
    !igy = nint(log(float(ny))/log(2.))
    !nxm = nx/2-1
    !nym = ny/2-1

   !$omp parallel do private(j, f, x) shared(zr, wfxy, wdy, nx, ny)
    do i=0,2*nx-1
        do j=0,ny-1
            f(j) = cmplx(zr(i,2*j),zr(i,2*j+1))
        enddo
        call ffttdo(f,ny,igy,wdy)
        call ftox(x,f,ny,wfxy)
        do j=0,nym
            if(mod(j,2).ne.0)then
                zr(i,j) = real(x(j))
            else
                zr(i,j) = imag(x(j))
            endif
        enddo
    enddo
  !$omp end parallel do
    
    !$omp parallel do private(i, f, x) shared(zr, z, wdx, wfxx, nx, ny)
    do j=0,nym
        do i=0,nx-1
            f(i) = cmplx(zr(2*i,j),zr(2*i+1,j))
        enddo
        call ffttdo(f,nx,igx,wdx)
        call ftox (x,f,nx,wfxx)
        if(mod(j,2).ne.0)then
            do i=0,nxm
                z(i,j) = x(i)
            enddo
        else
            do i=0,nxm
                z(i,j) = cmplx(-imag(x(i)),real(x(i)))
            enddo
        endif
    enddo
    !$omp end parallel do
    return
    end

    
    
! **********************************************************************
! subroutines specific for post processing
! **********************************************************************

    subroutine prefft_post(wo,n)
    use mod_main
    implicit none
    
    integer,intent(in) :: n
    complex*16,intent(out) :: wo(n/2-1)
    integer*2 :: irv
    integer :: i, ig, ipo
    real*8 :: argo
      
    ig=nint(log(float(n))/log2_v)
  
    do i=1,n/2-1
        ipo=irv(i,ig-1)
        argo=2.d0*pi*ipo/n
        wo(i)=cmplx(cos(argo),sin(argo))
    enddo
  
    return
    end
  
    ! **********************************************************************
  
    subroutine prefx_post(wy,n)
    use mod_main
    implicit none

    integer,intent(in) :: n
    complex*16,intent(out) :: wy(n/2-1)
    integer :: i
    real*8 :: arg

    do i=1,n/2-1
        arg=i*pi/n
        wy(i)=cmplx(-sin(arg),cos(arg))
    enddo
  
    return
    end
    