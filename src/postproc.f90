    subroutine postproc_map
    use mod_main
    implicit none

    integer :: i, j, nnx, nny, nnxm, nnym
    real*8 :: dx, dy, x, y
    real*8 :: pp1,pp2
  
    post_map = .true. !use to call the subroutine ORD
    
    !c0=6-2.5*log(2.5*ds)
    !c0=1./(c0*c0)
    !fs2=teta/c0
  
    nnx = id
    nny = jd
    nnxm = nnx/2-1
    nnym = nny/2-1
    do i=0,nnxm
        do j=0,nnym
            u(i,j) = (zero,zero)
            v(i,j) = (zero,zero)
            d(i,j) = (zero,zero)
            e(i,j) = (zero,zero)
        enddo
    enddo
  
    nxm = nx/2-1
    nym = ny/2-1
    
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_asc)
    open(50,status='old',file=trim(str_tmp))
    read(50,*) !t
    do i=0,nxm
        do j=0,nym
            read(50,*) pp1,pp2
            u(i,j) = cmplx(pp1,pp2)
            read(50,*) pp1,pp2
            v(i,j) = cmplx(pp1,pp2)
            read(50,*) pp1,pp2
            d(i,j) = cmplx(pp1,pp2)
            read(50,*) pp1,pp2
            e(i,j) = cmplx(pp1,pp2)
        enddo
    enddo
    close (50)
  
    ! prepare the coefficients for IFFT/FFT
    call prefft_post(wox,nnx)
    call prefft_post(woy,nny)
    call prefx_post(wxfx,nnx)
    call prefx_post(wxfy,nny)
    ! compute the IFFT    
    call atruod(u,ur,wxfx,wxfy,wox,woy,nnx,nny)
    call atrvod(v,vr,wxfx,wxfy,wox,woy,nnx,nny)
    call atruod(d,dr,wxfx,wxfy,wox,woy,nnx,nny)
    call atruod(e,er,wxfx,wxfy,wox,woy,nnx,nny)
    
        
    !     do i=0,2*nnx-1
    !         do j=nnym+nny+1,nnym+1,-1
    !             c=6+2.5*log(dr(i,j)/(2.5*ds))
    !             c=1./(c*c)
    !             u2v2s=sqrt(ur(i,j)*ur(i,j)+vr(i,j)*vr(i,j))
    !             tx=c*ur(i,j)*u2v2s/dr(i,j)
    !             ty=c*vr(i,j)*u2v2s/dr(i,j)
    !             teta=fs2*c*u2v2s*u2v2s
    !             er(i,j)=dr(i,j)
    !         enddo
    !     enddo

    dx = 2.d0*pi/(rlx*2*nnx)
    dy = 2.d0/nny
    
    ! maps of bottom elevation, velocity u&v, depth
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_u)
    open(61,file=trim(str_tmp),status='unknown')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_v)
    open(62,file=trim(str_tmp),status='unknown')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_d)
    open(63,file=trim(str_tmp),status='unknown')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_e)
    open(64,file=trim(str_tmp),status='unknown')
    do i=0,2*nnx-1
        x = i*dx
        write(61,*) '%i=',i
        write(62,*) '%i=',i
        write(63,*) '%i=',i
        write(64,*) '%i=',i
        do j=0,nny
            y = -1.d0+j*dy
            !MT write(60,*)x,y,er(i,j+nnym)
            write(61,*)x,y,ur(i,j+nny/2)
            write(62,*)x,y,vr(i,j+nny/2)
            write(63,*)x,y,dr(i,j+nny/2)
            write(64,*)x,y,er(i,j+nny/2)
        enddo
    enddo
    close(61)
    close(62)
    close(63)
    close(64)


    return
    end
    
    !********************************************************************
    !MT (Marco Toffolon): this is not needed any more because everything is done in Matlab
    !MT: I keep it as a memory of the previous procedure

    subroutine postproc_evol
    use mod_main
    implicit none

    integer :: i, j, nnx, nny, nnxm, nnym
    integer :: ii, jj
    real*8 :: fase
    real*8 :: pp1,pp2
  
    nnx=id
    nny=jd
    nnxm=nnx/2-1
    nnym=nny/2-1
    do i=0,nnxm
        do j=0,nnym
            u(i,j) = (zero,zero)
            v(i,j) = (zero,zero)
            d(i,j) = (zero,zero)
            e(i,j) = (zero,zero)
        enddo
    enddo
  
    nxm = nx/2-1
    nym = ny/2-1

    ii = evol_i
    jj = evol_j
    
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_interm_res)
    open(50,file=trim(str_tmp),status='old')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_interm_time)
    open(51,file=trim(str_tmp),status='old')
    
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_evol)
    open(60,file=trim(str_tmp),status='unknown')

    do while(.true.)
        read(51,*,err=100,end=101) t
        read(50,*,err=100,end=100) !t
        do i=0,nxm
            do j=0,nym
                read(50,*) pp1,pp2
                e(i,j) = cmplx(pp1,pp2)
            enddo
        enddo
        fase = atan2(imag(e(ii,jj)),real(e(ii,jj)))
        write(60,*) t,abs(e(ii,jj)),fase
    enddo

100 close (50)
101 close (51)    
    close (60)
  
    return
    end
    
