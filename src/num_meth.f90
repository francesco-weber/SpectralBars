! SUBROUTINE FOR NUMERICAL INTEGRATION
    
    subroutine solver
    use mod_main
    implicit none
    
    integer :: i, j
    complex*16 :: a13,a23,h,f1,f2,dudt,dvdt,dddt,dedt
    real*8 :: ud00,t1,t2
    complex*16 :: zt
    call cpu_time(t1)
!!$OMP PARALLEL 
    do i=0,nxm
        do j=0,nym
            a13 = 0.5d0*zilx(i)/f02
            a23 = 0.5d0*zily(j)/f02
            h = d(i,j)+e(i,j)
            f1 = 1.5d0*ax(i,j)-0.5d0*ax0(i,j)-a13*h
            f2 = 1.5d0*ay(i,j)-0.5d0*ay0(i,j)-a23*h
  
            dddt = -(zilx(i)*(1.5d0*ud(i,j)-0.5d0*ud0(i,j)) + zily(j)*(1.5d0*vd(i,j)-0.5d0*vd0(i,j)))
            dedt = -q0*(zilx(i)*(1.5d0*qx(i,j)-0.5d0*qx0(i,j)) + zily(j)*(1.5d0*qy(i,j)-0.5d0*qy0(i,j)))
            e(i,j) = e(i,j)+ dedt*dt
            d(i,j) = d(i,j)+ dddt*dt
            h = d(i,j)+e(i,j)
            dudt = f1-a13*h
            dvdt = f2-a23*h
            u(i,j) = u(i,j)+dudt*dt
            v(i,j) = v(i,j)+dvdt*dt
        enddo
    enddo
!!$omp end parallel 
    call cpu_time(t2)
    time_openmp = time_openmp + t2-t1
    !write(6,*) 'time = ', t2-t1

    u(0,0)= u(0,0)+beta*c0*dt
    v(0,0) = (zero,zero)
    e(0,0) = (zero,zero)
    ud00 = 0.d0
    do i=0,nxm
        do j=0,nym
            zt = u(i,j)*conjg(d(i,j))
            ud00 = ud00+real(zt)
        enddo
    enddo
    ud00 = ud00-real(u(0,0)*conjg(d(0,0)))
    d(0,0) = (1.d0-4.d0*ud00)/u(0,0)
  
    return
    end
    
    subroutine nl
        use mod_main
        implicit none
    
        integer :: i, j
        complex*16 :: zilyj
        real*8 :: c, u2v2, sqrtu2v2, tx, ty, teta, sind, cosd, dteta, fi, logdr
    
        complex*16 :: dudx(0:idm,0:jdm), dudy(0:idm,0:jdm), dvdx(0:idm,0:jdm), dvdy(0:idm,0:jdm), dedy(0:idm,0:jdm), dedx(0:idm,0:jdm)
        real*8 :: dudxr(0:id2,0:jd2), dudyr(0:id2,0:jd2), dvdxr(0:id2,0:jd2), dvdyr(0:id2,0:jd2), dedyr(0:id2,0:jd2), dedxr(0:id2,0:jd2), dedtaur(0:id2,0:jd2), dednir(0:id2,0:jd2)
        real*8 :: udr(0:id2,0:jd2), vdr(0:id2,0:jd2)
        real*8 :: axr(0:id2,0:jd2), ayr(0:id2,0:jd2), qxr(0:id2,0:jd2), qyr(0:id2,0:jd2)
    
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
    
        call atruod(u, ur, wxfx, wxfy, wox, woy, nx, ny)
        call atrvod(v, vr, wxfx, wxfy, wox, woy, nx, ny)
        call atruod(d, dr, wxfx, wxfy, wox, woy, nx, ny)
    
        ! Parallelizing the loop to check dr and process udr and vdr
        !$omp parallel do private(j) schedule(static)
        do j=0,2*ny-1
            do i=0,2*nx-1
                if(dr(i,j) < 2*ds) then
                    if(.not.prepare_stop) then
                        ! only the first time
                        write(6,*) 'Warning: depth < 0'
                        write(11,*) 'Warning: depth < 0 at t=',t
                        call out_error
                    endif
                    write(11,*) 'Warning: depth ', dr(i,j), ' < 0 at (i,j): ', i, j
                    prepare_stop = .true.
                endif
                udr(i,j) = ur(i,j) * dr(i,j)
                vdr(i,j) = vr(i,j) * dr(i,j)
            enddo
        enddo
        !$omp end parallel do
    
        call traudo(udr, ud, wfxx, wfxy, wdx, wdy, nx, ny)
        call travdo(vdr, vd, wfxx, wfxy, wdx, wdy, nx, ny)
    
        ! Combining dudx and dudy loops for better cache efficiency
        !$omp parallel do private(j) schedule(static)
        do j=0,nym
            do i=0,nxm
                dudx(i,j) = zilx(i) * u(i,j)
                dudy(i,j) = zily(j) * u(i,j)
            enddo
        enddo
        !$omp end parallel do
    
        call atruod(dudx, dudxr, wxfx, wxfy, wox, woy, nx, ny)
        call atrvod(dudy, dudyr, wxfx, wxfy, wox, woy, nx, ny)
    
        ! Parallelizing the loop to calculate axr
        !$omp parallel do private(j) schedule(static)
        do j=0,2*ny-1
            do i=0,2*nx-1
                axr(i,j) = -ur(i,j) * dudxr(i,j) - vr(i,j) * dudyr(i,j)
            enddo
        enddo
        !$omp end parallel do
    
        ! Combining dvdx, dvdy, dedy, and dedx loops for better cache efficiency
        !$omp parallel do private(j, zilyj) schedule(static)
        do j=0,nym
            zilyj = zily(j)
            do i=0,nxm
                dvdx(i,j) = zilx(i) * v(i,j)
                dvdy(i,j) = zilyj * v(i,j)
                dedy(i,j) = zilyj * e(i,j)
                dedx(i,j) = zilx(i) * e(i,j)
            enddo
        enddo
        !$omp end parallel do
    
        call atrvod(dvdx, dvdxr, wxfx, wxfy, wox, woy, nx, ny)
        call atruod(dvdy, dvdyr, wxfx, wxfy, wox, woy, nx, ny)
    
        ! Parallelizing the loop to calculate ayr
        !$omp parallel do private(j) schedule(static)
        do j=0,2*ny-1
            do i=0,2*nx-1
                ayr(i,j) = -ur(i,j) * dvdxr(i,j) - vr(i,j) * dvdyr(i,j)
            enddo
        enddo
        !$omp end parallel do
    
        call atrvod(dedy, dedyr, wxfx, wxfy, wox, woy, nx, ny)
        call atruod(dedx, dedxr, wxfx, wxfy, wox, woy, nx, ny)
    
        ! Minimizing redundant computations in the final loop
        !!$omp parallel do private(j, c, u2v2, sqrtu2v2, tx, ty, teta, sind, cosd, dteta, fi, logdr) shared(axr, ayr, qxr, qyr, ur, vr, dr, dedxr, dedyr)
        !$omp parallel do private(j) schedule(static)
        do j=0,2*ny-1
            do i=0,2*nx-1
                logdr = log(dr(i,j))
                c = 1.d0 / (6.d0 + 2.5d0 * logdr - dslog)**(2.d0)
                u2v2 = ur(i,j)*ur(i,j) + vr(i,j)*vr(i,j)
                sqrtu2v2 = sqrt(u2v2)
                tx = c * ur(i,j) * sqrtu2v2
                ty = c * vr(i,j) * sqrtu2v2
                axr(i,j) = axr(i,j) - beta * tx / dr(i,j)
                ayr(i,j) = ayr(i,j) - beta * ty / dr(i,j)
                teta = fs2 * c * u2v2
    
                sind = vr(i,j) / sqrtu2v2 - rsb / sqrt(teta) * ((ur(i,j) * dedyr(i,j) - vr(i,j) * dedxr(i,j)) / sqrtu2v2)
                cosd = sqrt(1.d0 - sind * sind)
    
                dteta = teta - (tetacr - 0.1 / beta * (beta * c0 * f02 - (ur(i,j) * dedxr(i,j) + vr(i,j) * dedyr(i,j))) / sqrtu2v2)
                if(dteta > 0.d0) then
                    fi = 8.d0 * dteta * sqrt(dteta)
                else
                    fi = 0.d0
                endif
                qxr(i,j) = fi * cosd
                qyr(i,j) = fi * sind
            enddo
        enddo
        !$omp end parallel do
    
        call traudo(axr, ax, wfxx, wfxy, wdx, wdy, nx, ny)
        call travdo(ayr, ay, wfxx, wfxy, wdx, wdy, nx, ny)
        call traudo(qxr, qx, wfxx, wfxy, wdx, wdy, nx, ny)
        call travdo(qyr, qy, wfxx, wfxy, wdx, wdy, nx, ny)
    
        return
    end
    

    !subroutine nl
    !    use mod_main
    !    implicit none
    !
    !    integer :: i, j
    !    complex*16 :: zilyj
    !    real*8 :: c, u2v2, sqrtu2v2, tx, ty, teta, sind, cosd, dteta, fi
    !
    !    complex*16 :: dudx(0:idm,0:jdm), dudy(0:idm,0:jdm), dvdx(0:idm,0:jdm), dvdy(0:idm,0:jdm), dedy(0:idm,0:jdm), dedx(0:idm,0:jdm)
    !    real*8 :: dudxr(0:id2,0:jd2), dudyr(0:id2,0:jd2), dvdxr(0:id2,0:jd2), dvdyr(0:id2,0:jd2), dedyr(0:id2,0:jd2), dedxr(0:id2,0:jd2), dedtaur(0:id2,0:jd2), dednir(0:id2,0:jd2)
    !    real*8 :: udr(0:id2,0:jd2), vdr(0:id2,0:jd2)
    !    real*8 :: axr(0:id2,0:jd2), ayr(0:id2,0:jd2), qxr(0:id2,0:jd2), qyr(0:id2,0:jd2)
    !
    !    ! Parallelizing the loop to copy data
    !    !$omp parallel do private(j) shared(ax0, ay0, ud0, vd0, qx0, qy0, ax, ay, ud, vd, qx, qy)
    !    do i=0,nxm
    !        do j=0,nym
    !            ax0(i,j) = ax(i,j)
    !            ay0(i,j) = ay(i,j)
    !            ud0(i,j) = ud(i,j)
    !            vd0(i,j) = vd(i,j)
    !            qx0(i,j) = qx(i,j)
    !            qy0(i,j) = qy(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call atruod(u, ur, wxfx, wxfy, wox, woy, nx, ny)
    !    call atrvod(v, vr, wxfx, wxfy, wox, woy, nx, ny)
    !    call atruod(d, dr, wxfx, wxfy, wox, woy, nx, ny)
    !
    !    ! Parallelizing the loop to check dr and process udr and vdr
    !    !$omp parallel do private(j) shared(udr, vdr, ur, vr, dr, prepare_stop, t)
    !    do i=0,2*nx-1
    !        do j=0,2*ny-1
    !            if(dr(i,j) < 2*ds) then
    !                if(.not.prepare_stop) then
    !                    ! only the first time
    !                    write(6,*) 'Warning: depth < 0'
    !                    write(11,*) 'Warning: depth < 0 at t=',t
    !                    call out_error
    !                endif
    !                write(11,*) 'Warning: depth ', dr(i,j), ' < 0 at (i,j): ', i, j
    !                prepare_stop = .true.
    !            endif
    !            udr(i,j) = ur(i,j) * dr(i,j)
    !            vdr(i,j) = vr(i,j) * dr(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call traudo(udr, ud, wfxx, wfxy, wdx, wdy, nx, ny)
    !    call travdo(vdr, vd, wfxx, wfxy, wdx, wdy, nx, ny)
    !
    !    ! Parallelizing the loop to calculate dudx and dudy
    !    !$omp parallel do private(j) shared(dudx, dudy, zilx, zily, u)
    !    do i=0,nxm
    !        do j=0,nym
    !            dudx(i,j) = zilx(i) * u(i,j)
    !            dudy(i,j) = zily(j) * u(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call atruod(dudx, dudxr, wxfx, wxfy, wox, woy, nx, ny)
    !    call atrvod(dudy, dudyr, wxfx, wxfy, wox, woy, nx, ny)
    !
    !    ! Parallelizing the loop to calculate axr
    !    !$omp parallel do private(j) shared(axr, ur, vr, dudxr, dudyr)
    !    do i=0,2*nx-1
    !        do j=0,2*ny-1
    !            axr(i,j) = -ur(i,j) * dudxr(i,j) - vr(i,j) * dudyr(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    ! Parallelizing the loop to calculate dvdx, dvdy, dedx, dedy
    !    !$omp parallel do private(j, zilyj) shared(dvdx, dvdy, dedy, dedx, zilx, zily, v, e)
    !    do j=0,nym
    !        zilyj = zily(j)
    !        do i=0,nxm
    !            dvdx(i,j) = zilx(i) * v(i,j)
    !            dvdy(i,j) = zilyj * v(i,j)
    !            dedy(i,j) = zilyj * e(i,j)
    !            dedx(i,j) = zilx(i) * e(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call atrvod(dvdx, dvdxr, wxfx, wxfy, wox, woy, nx, ny)
    !    call atruod(dvdy, dvdyr, wxfx, wxfy, wox, woy, nx, ny)
    !
    !    ! Parallelizing the loop to calculate ayr
    !    !$omp parallel do private(j) shared(ayr, ur, vr, dvdxr, dvdyr)
    !    do i=0,2*nx-1
    !        do j=0,2*ny-1
    !            ayr(i,j) = -ur(i,j) * dvdxr(i,j) - vr(i,j) * dvdyr(i,j)
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call atrvod(dedy, dedyr, wxfx, wxfy, wox, woy, nx, ny)
    !    call atruod(dedx, dedxr, wxfx, wxfy, wox, woy, nx, ny)
    !
    !    ! Parallelizing the final loop to calculate axr, ayr, qxr, qyr
    !    !$omp parallel do private(j, c, u2v2, sqrtu2v2, tx, ty, teta, sind, cosd, dteta, fi) shared(axr, ayr, qxr, qyr, ur, vr, dr, dedxr, dedyr)
    !    do i=0,2*nx-1
    !        do j=0,2*ny-1
    !            c = 1.d0/(6.d0 + 2.5d0 * log(dr(i,j)) - dslog)**(2.d0)
    !            u2v2 = ur(i,j)**2 + vr(i,j)**2
    !            sqrtu2v2 = sqrt(u2v2)
    !            tx = c * ur(i,j) * sqrtu2v2
    !            ty = c * vr(i,j) * sqrtu2v2
    !            axr(i,j) = axr(i,j) - beta * tx / dr(i,j)
    !            ayr(i,j) = ayr(i,j) - beta * ty / dr(i,j)
    !            teta = fs2 * c * u2v2
    !
    !            sind = vr(i,j) / sqrtu2v2 - rsb / sqrt(teta) * ((ur(i,j) * dedyr(i,j) - vr(i,j) * dedxr(i,j)) / sqrtu2v2)
    !            cosd = sqrt(1.d0 - sind * sind)
    !
    !            dteta = teta - (tetacr - 0.1 / beta * (beta * c0 * f02 - (ur(i,j) * dedxr(i,j) + vr(i,j) * dedyr(i,j))) / sqrtu2v2)
    !            if(dteta > 0.d0) then
    !                fi = 8.d0 * dteta * sqrt(dteta)
    !            else
    !                fi = 0.d0
    !            endif
    !            qxr(i,j) = fi * cosd
    !            qyr(i,j) = fi * sind
    !        enddo
    !    enddo
    !    !$omp end parallel do
    !
    !    call traudo(axr, ax, wfxx, wfxy, wdx, wdy, nx, ny)
    !    call travdo(ayr, ay, wfxx, wfxy, wdx, wdy, nx, ny)
    !    call traudo(qxr, qx, wfxx, wfxy, wdx, wdy, nx, ny)
    !    call travdo(qyr, qy, wfxx, wfxy, wdx, wdy, nx, ny)
    !
    !    return
    !end
    !
!
!********************************************************************

    !subroutine nl
    !use mod_main
    !implicit none
    !
    !integer :: i, j
    !complex*16 :: zilyj
    !real*8 :: c, u2v2, sqrtu2v2, tx, ty, teta, sind, cosd, dteta, fi
  !
    !complex*16 :: dudx(0:idm,0:jdm),dudy(0:idm,0:jdm),dvdx(0:idm,0:jdm),dvdy(0:idm,0:jdm),dedy(0:idm,0:jdm),dedx(0:idm,0:jdm)
    !real*8 :: dudxr(0:id2,0:jd2),dudyr(0:id2,0:jd2),dvdxr(0:id2,0:jd2),dvdyr(0:id2,0:jd2),dedyr(0:id2,0:jd2),dedxr(0:id2,0:jd2),dedtaur(0:id2,0:jd2),dednir(0:id2,0:jd2)
    !real*8 :: udr(0:id2,0:jd2),vdr(0:id2,0:jd2)
    !real*8 :: axr(0:id2,0:jd2),ayr(0:id2,0:jd2),qxr(0:id2,0:jd2),qyr(0:id2,0:jd2)
!
    !
    !!equivalence (dudx,dvdx),(dudxr,dvdxr,qxr,udr)
    !!equivalence (dudy,dvdy),(dudyr,dvdyr,qyr,vdr)
  !
    !do i=0,nxm
    !    do j=0,nym
    !        ax0(i,j) = ax(i,j)
    !        ay0(i,j) = ay(i,j)
    !        ud0(i,j) = ud(i,j)
    !        vd0(i,j) = vd(i,j)
    !        qx0(i,j) = qx(i,j)
    !        qy0(i,j) = qy(i,j)
    !    enddo
    !enddo
  !
    !call atruod(u,ur,wxfx,wxfy,wox,woy,nx,ny)
    !call atrvod(v,vr,wxfx,wxfy,wox,woy,nx,ny)
    !call atruod(d,dr,wxfx,wxfy,wox,woy,nx,ny)
    !do i=0,2*nx-1
    !    do j=0,2*ny-1
    !        if(dr(i,j).lt.2*ds) then
    !            if(.not.prepare_stop) then
    !                ! only the first time
    !                write(6,*) 'Warning: depth < 0'
    !                write(11,*) 'Warning: depth < 0 at t=',t
    !                call out_error
    !            endif
    !            write(11,*) 'Warning: depth ',dr(i,j),' < 0 at (i,j): ',i,j
    !            prepare_stop = .true.
    !        endif
    !        udr(i,j) = ur(i,j)*dr(i,j)
    !        vdr(i,j) = vr(i,j)*dr(i,j)
    !    enddo
    !enddo
    !call traudo(udr,ud,wfxx,wfxy,wdx,wdy,nx,ny)
    !call travdo(vdr,vd,wfxx,wfxy,wdx,wdy,nx,ny)
    !do i=0,nxm
    !    do j=0,nym
    !        dudx(i,j) = zilx(i)*u(i,j)
    !        dudy(i,j) = zily(j)*u(i,j)
    !    enddo
    !enddo
    !call atruod(dudx,dudxr,wxfx,wxfy,wox,woy,nx,ny)
    !call atrvod(dudy,dudyr,wxfx,wxfy,wox,woy,nx,ny)
    !do i=0,2*nx-1
    !    do j=0,2*ny-1
    !        axr(i,j) = -ur(i,j)*dudxr(i,j)-vr(i,j)*dudyr(i,j)
    !    enddo
    !enddo
    !do j=0,nym
    !    zilyj = zily(j)
    !    do i=0,nxm
    !        dvdx(i,j) = zilx(i)*v(i,j)
    !        dvdy(i,j) = zilyj*v(i,j)
    !        dedy(i,j) = zilyj*e(i,j)
    !        dedx(i,j) = zilx(i)*e(i,j)
    !    enddo
    !enddo
    !call atrvod(dvdx,dvdxr,wxfx,wxfy,wox,woy,nx,ny)
    !call atruod(dvdy,dvdyr,wxfx,wxfy,wox,woy,nx,ny)
    !do i=0,2*nx-1
    !    do j=0,2*ny-1
    !        ayr(i,j) = -ur(i,j)*dvdxr(i,j)-vr(i,j)*dvdyr(i,j)
    !    enddo
    !enddo
!
    !call atrvod(dedy,dedyr,wxfx,wxfy,wox,woy,nx,ny)
    !call atruod(dedx,dedxr,wxfx,wxfy,wox,woy,nx,ny)
!
    !do i=0,2*nx-1
    !    do j=0,2*ny-1
    !        c = 1.d0/(6.d0+2.5d0*log(dr(i,j))-dslog)**(2.d0)
    !        u2v2 = ur(i,j)*ur(i,j)+vr(i,j)*vr(i,j)
    !        sqrtu2v2=sqrt(u2v2)
    !        tx = c*ur(i,j)*sqrtu2v2
    !        ty = c*vr(i,j)*sqrtu2v2
    !        axr(i,j) = axr(i,j)-beta*tx/dr(i,j)
    !        ayr(i,j) = ayr(i,j)-beta*ty/dr(i,j)
    !        teta = fs2*c*u2v2
!
    !        sind = vr(i,j)/sqrtu2v2-rsb/sqrt(teta)*((ur(i,j)*dedyr(i,j)-vr(i,j)*dedxr(i,j))/sqrtu2v2)
    !        cosd = sqrt(1.d0-sind*sind)
!
    !        dteta = teta-(tetacr-0.1/(beta)*(beta*c0*f02-(ur(i,j)*dedxr(i,j)+vr(i,j)*dedyr(i,j)))/sqrtu2v2)
    !        if(dteta.gt.0.d0) then
    !            fi = 8.d0*(dteta)*sqrt(dteta)
    !        else
    !            fi = 0.d0
    !        endif
    !        qxr(i,j) = fi*cosd
    !        qyr(i,j) = fi*sind
    !    enddo
    !enddo
    !call traudo(axr,ax,wfxx,wfxy,wdx,wdy,nx,ny)
    !call travdo(ayr,ay,wfxx,wfxy,wdx,wdy,nx,ny)
    !call traudo(qxr,qx,wfxx,wfxy,wdx,wdy,nx,ny)
    !call travdo(qyr,qy,wfxx,wfxy,wdx,wdy,nx,ny)
  !
    !return
    !end
    !
!