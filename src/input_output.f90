    subroutine input
    use mod_main
    implicit none
    integer :: i, j
    real(8) :: tmax, tout
      
    open(40,status='old',file='input.txt') !,err=100)
    read(40,*) !% this is the input file
    read(40,*) !% physical parameters
    read(40,*) teta0
    read(40,*) ds
    read(40,*) beta
    read(40,*) rlx
    read(40,*) rt
    read(40,*) !% numerical parameters
    read(40,*) nx
    read(40,*) ny
    read(40,*) dt
    read(40,*) !% initial random perturbation
    read(40,*) eps
    read(40,*) random_seed
    read(40,*) !% parameters for the output
    read(40,*) tmax
    read(40,*) tout
    close(40)
    
    open(40,status='old',file='input_filenames.txt') !,err=100)
    read(40,*) !% this is the input file
    read(40,*) working_folder
    read(40,*) !% file names
    read(40,*) file_res_bin
    read(40,*) file_res_asc
    read(40,*) file_interm_res
    read(40,*) file_interm_time
    read(40,*) !% post processing
    read(40,*) file_res_map_u
    read(40,*) file_res_map_v
    read(40,*) file_res_map_d
    read(40,*) file_res_map_e
    read(40,*) file_res_evol
    read(40,*) evol_i
    read(40,*) evol_j
    close(40)
    
    ntmax = nint(tmax/dt)
    nusc  = nint(tout/dt)
    !log2_v = log(2.)
    !ig  = nint(log(float(n))/log(2.d0))
    igx = nint(log(float(nx))/log2_v)
    igy = nint(log(float(ny))/log2_v)
    nxm = nx/2-1
    nym = ny/2-1
    
    return
100  pause 'Error on input files'
    stop
    end
    
    !********************************************************************
    
    subroutine initial_conditions
    use ifport !mt
    use mod_main
    implicit none
    integer :: i, j
    
    ! set random initial conditions
    call seed(random_seed) ! initialize the random seed

    t0 = 0
    do i=0,nxm
        do j=0,nym
            !mt
            if(i.eq.0.and.j.eq.0) then
                u(i,j) = 1.d0
                d(i,j) = 1.d0
            else
                u(i,j) = 0
                d(i,j) = 0
            endif
            
            !MT start of modified code: the first harmonic is fixed (!!!)
            if(i.eq.1.and.j.eq.1.and.eps.eq.0.d0) then
                e(i,j) = 0.01d0    ! fixed initial amplitude
            else
                if(i.lt.20.and.j.lt.20) then
                    e(i,j) = random(0)*eps
                    !e(i,j) = eps
                    !write(*,*) 'ATTIVARE RANDOM SEED'
                else
                     e(i,j) = 0
                endif
            endif
            !MT end of modified code
            
            v(i,j) = 0
            ax(i,j) = 0
            ay(i,j) = 0
            ud(i,j) = 0
            vd(i,j) = 0
            qx(i,j) = 0
            qy(i,j) = 0
            ax0(i,j) = 0
            ay0(i,j) = 0
            ud0(i,j) = 0
            vd0(i,j) = 0
            qx0(i,j) = 0
            qy0(i,j) = 0
        enddo
    enddo
      
    return
    end
    
    !********************************************************************
  
    subroutine out_intermediate
    use mod_main, only: t, nxm, nym, u, v, d, e, ureal, vreal, dreal, ereal, wox,woy,id,jd, wxfx, wxfy, wdx,wdy ,nx,ny,wfxx,wfxy,pi,rlx, str_tmp, file_res_map_u, file_res_map_v, file_res_map_d, file_res_map_e, working_folder,post_map
    implicit none
    
    integer :: i, j, nnx, nny
    real*8 :: dx, dy, x, y

    post_map= .true.



      
    !  scrittura su file risultati intermedi
  
    !fase=atan2(imag(e(1,1)),real(e(1,1)))
    !write(60,*) t,abs(e(1,1)),fase

    write(60,*) '%time=', t
    do i=0,nxm
        do j=0,nym
         !write(60,*) u(i,j)
         !write(60,*) v(i,j)
         !write(60,*) d(i,j)
            write(60,*) real(e(i,j)),imag(e(i,j)),i,j
        enddo
    enddo
    !write 
    write(61,*) t

    !post_map = .false. !use to call the subroutine ORD



    nnx = id
    nny = jd


    ! prepare the coefficients for IFFT/FFT
    call prefft_post(wox,nnx)
    call prefft_post(woy,nny)
    call prefx_post(wxfx,nnx)
    call prefx_post(wxfy,nny)
    !! compute the IFFT    
    call atruod(u,ureal,wxfx,wxfy,wox,woy,nnx,nny)
    call atrvod(v,vreal,wxfx,wxfy,wox,woy,nnx,nny)
    call atruod(d,dreal,wxfx,wxfy,wox,woy,nnx,nny)
    call atruod(e,ereal,wxfx,wxfy,wox,woy,nnx,nny)
    
    !call prefft(wdx,wox,nx)
    !call prefft(wdy,woy,ny)
    !call prefx(wfxx,wxfx,nx)
    !call prefx(wfxy,wxfy,ny)

    dx = 2.d0*pi/(rlx*2*nnx)
    dy = 2.d0/nny
    

    ! maps of bottom elevation, velocity u&v, depth
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_u)
    open(62,file=trim(str_tmp),status='unknown',position='append')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_v)
    open(63,file=trim(str_tmp),status='unknown',position='append')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_d)
    open(64,file=trim(str_tmp),status='unknown',position='append')
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_e)
    open(65,file=trim(str_tmp),status='unknown',position='append')

    write(62,*)'%time=', t
    write(63,*)'%time=', t
    write(64,*)'%time=', t
    write(65,*)'%time=', t

    do i=0,2*nnx-1
        x = i*dx
        write(62,*) '%i=',i
        write(63,*) '%i=',i
        write(64,*) '%i=',i
        write(65,*) '%i=',i
        do j=0,nny
            y = -1.d0+j*dy
            write(62,*)x,y,ureal(i,j+nny/2)
            write(63,*)x,y,vreal(i,j+nny/2)
            write(64,*)x,y,dreal(i,j+nny/2)
            write(65,*)x,y,ereal(i,j+nny/2)
        enddo
    enddo
    close(62)
    close(63)
    close(64)
    close(65)

    post_map = .false. !use to call the subroutine ORD


    return
    end
  
    !********************************************************************
  
    subroutine out_final
    use mod_main
    implicit none
    
    integer :: i, j
  
    !  scrittura su file risultati finali
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_bin)
    open(70,file=trim(str_tmp),form='unformatted')
    write(70) t
    do i=0,nxm
        do j=0,nym
            write(70) u(i,j)
            write(70) v(i,j)
            write(70) d(i,j)
            write(70) e(i,j)
            write(70) ax(i,j)
            write(70) ay(i,j)
            write(70) ud(i,j)
            write(70) vd(i,j)
            write(70) qx(i,j)
            write(70) qy(i,j)
            write(70) ax0(i,j)
            write(70) ay0(i,j)
            write(70) ud0(i,j)
            write(70) vd0(i,j)
            write(70) qx0(i,j)
            write(70) qy0(i,j)
        enddo
    enddo
    close(70)

    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_asc)
    open(80,file=trim(str_tmp))
    write(80,*) '%time=', t
    do i=0,nxm
        do j=0,nym
            write(80,*) real(u(i,j)),imag(u(i,j)),i,j
            write(80,*) real(v(i,j)),imag(v(i,j)),i,j
            write(80,*) real(d(i,j)),imag(d(i,j)),i,j
            write(80,*) real(e(i,j)),imag(e(i,j)),i,j
        enddo
    enddo
    close(80)
  
    return
    end
    
    
    
    
    !********************************************************************
    ! write the files when an error occurs
    subroutine out_error
    use mod_main
    implicit none
    integer :: i, j
  
    ! write the amplitudes
    write(str_tmp,'(3a)') trim(working_folder),'/','err_nan.txt'
    open(80,file=trim(str_tmp),status='unknown')
    write(80,*) '% amplitudes (u,v,d,e) at time=', t
    do i=0,nxm
        do j=0,nym
            write(80,*) real(u(i,j)),imag(u(i,j)),i,j
            write(80,*) real(v(i,j)),imag(v(i,j)),i,j
            write(80,*) real(d(i,j)),imag(d(i,j)),i,j
            write(80,*) real(e(i,j)),imag(e(i,j)),i,j
        enddo
    enddo
    close(80)
    
    ! write the real values before the Nan occurred
    write(str_tmp,'(3a)') trim(working_folder),'/','err_map_dr.txt'
    open(80,file=trim(str_tmp),status='unknown')
    write(80,*) '% un-ordered map (d) at time=', t
    do i=0,2*nx-1
        do j=0,2*ny-1
            write(80,*) dr(i,j),i,j
        enddo
    enddo
    close(80)
    
    !! write the real values before the Nan occurred
    !write(str_tmp,'(3a)') trim(working_folder),'/','err_map.txt'
    !open(80,file=trim(str_tmp),status='unknown')
    !write(80,*) '% un-ordered map (e, d) at time=', t
    !do i=0,2*nx-1
    !    do j=0,2*ny-1
    !        write(80,*) er_tmp(i,j),dr_tmp(i,j),i,j
    !    enddo
    !enddo
    !close(80)
    
    return
    end
