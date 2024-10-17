    program main
    use mod_main
    USE OMP_LIB
    implicit none
    
    integer*4 :: it
    integer :: i, j
    real*8 :: maxmode
    real*8 :: t_tot_start, t_tot_end
    integer :: imaxmode, jmaxmode
    
    prepare_stop = .false. ! to stop the execution after an error


    time_openmp = 0.d0
    call input
    call allocate_vars
    call param
    call initial_conditions
    
    post_map = .false. ! avoid calling ORD in atruod and atrvod
    write(str_tmp,'(3a)') trim(working_folder),'/','ctrl.txt'
    open(11,file=trim(str_tmp),status='unknown')
    write(11,*) 'start'
    write(11,'(2a)') 'working folder: ',trim(working_folder)
    write(*,'(2a)') 'working folder: ',trim(working_folder)


    call prefft(wdx,wox,nx)
    call prefft(wdy,woy,ny)
    call prefx(wfxx,wxfx,nx)
    call prefx(wfxy,wxfy,ny)
  
    if(t0.eq.zero) then
        ! first step with Eulero
        call nl
        !$omp parallel do private(j) schedule(static)
        do j=0,nym
            do i=0,nxm
                ax0(i,j) = ax(i,j)
                ay0(i,j) = ay(i,j)
                ud0(i,j) = ud(i,j)
                vd0(i,j) = vd(i,j)
                qx0(i,j) = qx(i,j)
                qy0(i,j) = qy(i,j)
            enddo
        enddo
        !$omp end parallel do
        call solver
        call nl
    endif
 
    ! main loop on time
    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_interm_time)
    open(61,file=trim(str_tmp),status='unknown') ! open file for time of intermediate results

    write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_interm_res)
    open(60,file=trim(str_tmp),status='unknown') ! open file for intermediate results

    !write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_u)
    !open(62,file=trim(str_tmp),status='unknown',position='append')
    !write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_v)
    !open(63,file=trim(str_tmp),status='unknown',position='append')
    !write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_d)
    !open(64,file=trim(str_tmp),status='unknown',position='append')
    !write(str_tmp,'(3a)') trim(working_folder),'/',trim(file_res_map_e)
    !open(65,file=trim(str_tmp),status='unknown',position='append')




    t = t0
    !call out_intermediate   ! write initial conditions
    
    write(str_tmp,'(3a)') trim(working_folder),'/','res_max_ampl.txt'
    open(12,file=trim(str_tmp),status='unknown') ! values of the max amplitude
    write(12,*) '% max amplitude: t, A, i, j'

    call cpu_time(t_tot_start)

    do it=1,ntmax

        t = t0 + it*dt
        !call atruod(e,er_tmp,wxfx,wxfy,wox,woy,nx,ny) !mt   temp -- TO BE REMOVED
        !call atruod(d,dr_tmp,wxfx,wxfy,wox,woy,nx,ny) !mt   temp -- TO BE REMOVED
        call solver
        call nl
        if(prepare_stop) goto 999
        !mt
        if(u(1,1).ne.u(1,1)) then ! strange way to define a 'not a number'
            ! if there is a NaN, write the results obtained so far and stop the computation
            prepare_stop = .true.
            print*,'NaN detected: stop execution at it=',it
            write(11,*) 'NaN detected: stop execution at it=',it
            write(11,*) 'corresponing to t=',t
            call out_error
            ! recover the last saved previous state
            u = up
            v = vp
            d = dp
            e = ep
            t = tp
            goto 999
        else
            ! otherwise update the temporary variables
            up = u
            vp = v
            dp = d
            ep = e
            tp = t
        endif
        if(mod(it,nusc).eq.0) then
            write(*,*) 'out',it,'/',ntmax !mt
            write(11,*) 'out',it,'/',ntmax,' t=',t
            call out_intermediate
            ! find the mode with max amplitude
            maxmode = zero
            do i=0,nxm
                do j=0,nym
                    if(abs(e(i,j)).gt.maxmode) then
                        maxmode = abs(e(i,j))
                        imaxmode = i
                        jmaxmode = j
                    endif
                enddo
            enddo
            write(11,*) 'maxmode ',maxmode,' at i,j:',imaxmode,jmaxmode
            write(12,*) t,maxmode,imaxmode,jmaxmode
        endif
    enddo
    call cpu_time(t_tot_end)


999 continue
    close(60)
    close(61)
    close(12)

    write(6,*) 'time_openmp = ', time_openmp
    write(6,*) 'time_tot = ', t_tot_end - t_tot_start
    
    call out_final
    
    print*,'end of computation'
    print*,'-'
    print*,'postprocessing: map'
    !call postproc_map
    
    print*,'postprocessing: evol'
    print*,'harmonic:',evol_i,evol_j
    call postproc_evol
    
    print*,'end program'
    !pause   ! this can be removed
    
    write(11,*) 'end'
    close(11)
    
    stop
    end