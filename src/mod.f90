    module mod_main
    implicit none

    ! Variables
    complex*16,allocatable,dimension(:,:) :: u, v, d, e, ax, ax0, ay, ay0, ud, ud0, vd, vd0, qx, qx0, qy, qy0
    complex*16,allocatable,dimension(:) :: wxfx, wxfy, wfxx, wfxy, wox, woy, wdx, wdy
    complex*16,allocatable,dimension(:) :: zilx,zily
    complex*16 :: zi
    real*8 :: dslog, fs2, tetacr, f02, q0, rt, c0, rsb

    integer*4 :: ntmax
    integer*4 :: igx,igy

    real*8 :: t0, t, dt, teta0, ds, beta, rlx
    integer :: nx, ny, nxm, nym, ilu, nusc
    !integer :: igx, igy

    integer :: id, jd, ijd, idm, jdm, ijdm, id2, jd2
    
    real*8,parameter :: zero = 0.d0
    real*8,parameter :: log2_v = log(2.0)

    !mt
    complex*16,allocatable,dimension(:,:) :: up, vp, dp, ep
    real*8,allocatable,dimension(:,:) :: ur,vr,dr,er,ureal,vreal,dreal,ereal

    real*8 :: tp

    real*8 :: pi

    real*8 :: time_openmp

    real*8 :: eps
    integer :: random_seed
    
    character(128) :: working_folder
    character(256) :: str_tmp
    character(128) :: file_res_bin, file_res_asc, file_interm_res, file_interm_time, file_res_evol
    character(128) :: file_res_map_u, file_res_map_v, file_res_map_d, file_res_map_e
    integer :: evol_i, evol_j
    
    logical :: post_map, prepare_stop

    real*8,allocatable :: er_tmp(:,:),dr_tmp(:,:)

    end module
