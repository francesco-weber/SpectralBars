    subroutine allocate_vars
    use mod_main

    implicit none
    
    id = nx
    jd = ny
    ijd = max(id,jd)
    idm = id/2-1
    jdm = jd/2-1
    id2 = 2*id-1
    jd2 = 2*jd-1
    ijdm = ijd/2-1
    
    nxm=nx/2-1
    nym=ny/2-1

    
    allocate(u(0:idm,0:jdm),v(0:idm,0:jdm),d(0:idm,0:jdm),e(0:idm,0:jdm))
    allocate(ax(0:idm,0:jdm),ax0(0:idm,0:jdm),ay(0:idm,0:jdm),ay0(0:idm,0:jdm))
    allocate(ud(0:idm,0:jdm),ud0(0:idm,0:jdm),vd(0:idm,0:jdm),vd0(0:idm,0:jdm))
    allocate(qx(0:idm,0:jdm),qx0(0:idm,0:jdm),qy(0:idm,0:jdm),qy0(0:idm,0:jdm))
    allocate(wxfx(idm),wxfy(jdm),wfxx(idm),wfxy(jdm))
    allocate(wox(idm),woy(jdm),wdx(idm),wdy(jdm))
    allocate(zilx(0:idm),zily(0:jdm))
    
    ! temporary variables to store the results at a time before NaN
    allocate(up(0:idm,0:jdm),vp(0:idm,0:jdm),dp(0:idm,0:jdm),ep(0:idm,0:jdm))

    !allocate(er_tmp(0:id2,0:jd2),dr_tmp(0:id2,0:jd2))
    
    allocate(ur(0:id2,0:jd2),vr(0:id2,0:jd2),dr(0:id2,0:jd2),er(0:id2,0:jd2),ureal(0:id2,0:jd2),vreal(0:id2,0:jd2),dreal(0:id2,0:jd2),ereal(0:id2,0:jd2))

    return    
    end
    