module file_IO
    use iso_fortran_env, only: stdout => output_unit, dp => real64, sp => real32
    implicit none

    private
    public :: read_target_value, write_eam_file, generate_eam_parameter, summary, update_eam_parameter
    character(len=*), parameter, private :: filename = "../parameter.nml"
    integer, parameter, private :: nr = 2000    !< r划分多少份
    integer, parameter, private :: nrho = 2000  !< rho划分多少份
    real(dp), parameter, private :: pi = 4.0d0*atan(1.0d0)  !< 圆周率pi
    logical, private :: first_read_parameter = .False.
    ! 以下是势函数参数
    real(dp), private, save :: re 
    real(dp), private, save :: fe 
    real(dp), private, save :: rou_e 
    real(dp), private, save :: rou_s 
    real(dp), private, save :: alpha 
    real(dp), private, save :: beta 
    real(dp), private, save :: A 
    real(dp), private, save :: B 
    real(dp), private, save :: cai 
    real(dp), private, save :: lamda 
    real(dp), private, save :: Fn0 
    real(dp), private, save :: Fn1 
    real(dp), private, save :: Fn2 
    real(dp), private, save :: Fn3 
    real(dp), private, save :: F0 
    real(dp), private, save :: F1 
    real(dp), private, save :: F2 
    real(dp), private, save :: F3 
    real(dp), private, save :: eta 
    real(dp), private, save :: F_e 
    real(dp), private, save :: rhoin
    real(dp), private, save :: rhoout
    ! 以下是规定了势函数参数上下限的参数
    real(dp), private, dimension(2), save :: re_bound        ! re参数的取值上下限
    real(dp), private, dimension(2), save :: fe_bound        ! fe参数的取值上下限 
    real(dp), private, dimension(2), save :: rou_e_bound    ! rou_e参数的取值上下限
    real(dp), private, dimension(2), save :: rou_s_bound    ! rou_s参数的取值上下限
    real(dp), private, dimension(2), save :: alpha_bound   ! alpha参数的取值上下限
    real(dp), private, dimension(2), save :: beta_bound       ! beta参数的取值上下限
    real(dp), private, dimension(2), save :: A_bound    ! A参数的取值上下限
    real(dp), private, dimension(2), save :: B_bound         ! B参数的取值上下限
    real(dp), private, dimension(2), save :: cai_bound       ! cai参数的取值上下限
    real(dp), private, dimension(2), save :: lamda_bound     ! lamda参数的取值上下限
    real(dp), private, dimension(2), save :: Fn0_bound      ! Fn0参数的取值上下限
    real(dp), private, dimension(2), save :: Fn1_bound      ! Fn1参数的取值上下限
    real(dp), private, dimension(2), save :: Fn2_bound       ! Fn2参数的取值上下限
    real(dp), private, dimension(2), save :: Fn3_bound      ! Fn3参数的取值上下限
    real(dp), private, dimension(2), save :: F0_bound       ! F0参数的取值上下限
    real(dp), private, dimension(2), save :: F1_bound        ! F1参数的取值上下限
    real(dp), private, dimension(2), save :: F2_bound        ! F2参数的取值上下限
    real(dp), private, dimension(2), save :: F3_bound       ! F3参数的取值上下限
    real(dp), private, dimension(2), save :: eta_bound      ! eta参数的取值上下限
    real(dp), private, dimension(2), save :: F_e_bound     ! F_e参数的取值上下限
    ! 迭代过程中的临时势函数参数
    ! 以下是势函数参数
    real(dp), private, save :: re_ 
    real(dp), private, save :: fe_
    real(dp), private, save :: rou_e_ 
    real(dp), private, save :: rou_s_ 
    real(dp), private, save :: alpha_ 
    real(dp), private, save :: beta_ 
    real(dp), private, save :: A_ 
    real(dp), private, save :: B_ 
    real(dp), private, save :: cai_ 
    real(dp), private, save :: lamda_ 
    real(dp), private, save :: Fn0_ 
    real(dp), private, save :: Fn1_ 
    real(dp), private, save :: Fn2_ 
    real(dp), private, save :: Fn3_ 
    real(dp), private, save :: F0_ 
    real(dp), private, save :: F1_ 
    real(dp), private, save :: F2_ 
    real(dp), private, save :: F3_ 
    real(dp), private, save :: eta_ 
    real(dp), private, save :: F_e_
contains

    ! 此函数用于生成区间在[lo, hi]之间随机均匀分布的随机数
    function rand_normal(mu, sigma) result(res)
        real(dp), intent(in) :: mu, sigma
        real(dp) :: res
        real(dp) :: u1, u2
        real(dp) :: z0

        ! 生成2个0-1之间均匀分布的随机数
        call random_number(u1)
        call random_number(u2)
        ! 利用 Box - Muller 变换生成标准正态分布随机数
        z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
        ! 线性缩放, 生成均质为mu, 标准差为sigma的正态分布的随机数
        res = mu + sigma * z0
    end function rand_normal

    subroutine summary(cost)
        real(dp), intent(in) :: cost

        write(stdout, '(A)') "Iteration finished!"
        write(stdout, '(A, F24.16)') 're =', re 
        write(stdout, '(A, F24.16)') 'fe =', fe 
        write(stdout, '(A, F24.16)') 'rou_e =', rou_e 
        write(stdout, '(A, F24.16)') 'rou_s =', rou_s 
        write(stdout, '(A, F24.16)') 'alpha =', alpha
        write(stdout, '(A, F24.16)') 'beta =', beta
        write(stdout, '(A, F24.16)') 'A =', A
        write(stdout, '(A, F24.16)') 'B =', B 
        write(stdout, '(A, F24.16)') 'cai =', cai
        write(stdout, '(A, F24.16)') 'lamda =', lamda 
        write(stdout, '(A, F24.16)') 'Fn0 =', Fn0
        write(stdout, '(A, F24.16)') 'Fn1 =', Fn1 
        write(stdout, '(A, F24.16)') 'Fn2 =', Fn2
        write(stdout, '(A, F24.16)') 'Fn3 =', Fn3
        write(stdout, '(A, F24.16)') 'F0 =', F0 
        write(stdout, '(A, F24.16)') 'F1 =', F1
        write(stdout, '(A, F24.16)') 'F2 =', F2
        write(stdout, '(A, F24.16)') 'F3 =', F3 
        write(stdout, '(A, F24.16)') 'eta =', eta 
        write(stdout, '(A, F24.16)') 'F_e =', F_e
        write(stdout, '(A, es24.16)') 'The final cost is: ', cost
    end subroutine

    ! 在用户给定的定义域之内更新eam势函数参数
    subroutine generate_eam_parameter()
        Fn0 = Fn0 + rand_normal(0.0d0, 0.1d0)
        ! 约束Fn0在固定范围之内
        if (Fn0 > Fn0_bound(2)) Fn0 = Fn0_bound(2)
        if (Fn0 < Fn0_bound(1)) Fn0 = Fn0_bound(1)
        Fn1 = Fn1 + rand_normal(0.0d0, 0.1d0)
        ! 约束Fn1在固定范围之内
        if (Fn1 > Fn1_bound(2)) Fn1 = Fn1_bound(2)
        if (Fn1 < Fn1_bound(1)) Fn1 = Fn1_bound(1)
        Fn2 = Fn2 + rand_normal(0.0d0, 0.1d0)
        ! 约束Fn2在固定范围之内
        if (Fn2 > Fn2_bound(2)) Fn2 = Fn2_bound(2)
        if (Fn2 < Fn2_bound(1)) Fn2 = Fn2_bound(1)
        Fn3 = Fn3 + rand_normal(0.0d0, 0.1d0)
        ! 约束Fn3在固定范围之内
        if (Fn3 > Fn3_bound(2)) Fn3 = Fn3_bound(2)
        if (Fn3 < Fn3_bound(1)) Fn3 = Fn3_bound(1)
    end subroutine generate_eam_parameter

    subroutine update_eam_parameter(flag)
        integer, intent(in) :: flag
        if (flag == 1) then
            re      = re_
            fe      = fe_ 
            rou_e   = rou_e_
            rou_s   = rou_s_
            alpha   = alpha_
            beta    = beta_ 
            A       = A_
            B       = B_
            cai     = cai_
            lamda   = lamda_
            Fn0     = Fn0_
            Fn1     = Fn1_
            Fn2     = Fn2_
            Fn3     = Fn3_
            F0      = F0_
            F1      = F1_
            F2      = F2_
            F3      = F3_
            eta     = eta_
            F_e     = F_e_
        else 
            re_      = re
            fe_      = fe
            rou_e_   = rou_e
            rou_s_   = rou_s
            alpha_   = alpha
            beta_    = beta
            A_       = A
            B_       = B
            cai_     = cai
            lamda_   = lamda
            Fn0_     = Fn0
            Fn1_     = Fn1
            Fn2_     = Fn2
            Fn3_     = Fn3
            F0_      = F0
            F1_      = F1
            F2_      = F2
            F3_      = F3
            eta_     = eta
            F_e_     = F_e
        end if
    end subroutine update_eam_parameter

    ! 读取拟合目标, 以及拟合参数上下限
    subroutine read_target_value(lat, a0, Ec, Efv, C11, C12, C44, Bulk)
        character(len=3), intent(out) :: lat    !< 晶格类型
        real(dp), intent(out) :: a0             !< 晶格常数
        real(dp), intent(out) :: Ec             !< 内聚能
        real(dp), intent(out) :: Efv            !< 空位形成能
        real(dp), intent(out) :: C11            !< 弹性常数C11
        real(dp), intent(out) :: C12            !< 弹性常数C22
        real(dp), intent(out) :: C44            !< 弹性常数C44
        real(dp), intent(out) :: Bulk            !< 体积模量
        integer :: fileid, iosuccess
        logical :: file_exist
        namelist/targetvalue/ lat, a0, Ec, Efv, C11, C12, C44, Bulk
        namelist/bound/ re_bound, fe_bound, rou_e_bound, rou_s_bound, alpha_bound,  &
        beta_bound, A_bound, B_bound, cai_bound, lamda_bound, Fn0_bound, Fn1_bound, &
        Fn2_bound, Fn3_bound, F0_bound, F1_bound, F2_bound, F3_bound, eta_bound, F_e_bound

        inquire(file=filename, exist=file_exist)
        if (.not. file_exist) then
            write(stdout, "(A)") "File ../parameter.nml doesn't exist!"
            stop
        end if
        open(newunit=fileid, file=filename,action="read")
        read(fileid, nml=targetvalue, iostat=iosuccess)
        if (iosuccess /= 0) then
            write(stdout, "(A)") "Error occuring during read ../parameter.nml, please check your format."
            write(stdout, "(A)") "The reference format of ../parameter.nml like:"
            write(stdout, "(A)") '&targetvalue'
            write(stdout, "(A)") 'lat = "bcc"'
            write(stdout, "(A)") 'a0 = 3.03D0'
            write(stdout, "(A)") 'Ec = 5.31D0'
            write(stdout, "(A)") 'Efv = 2.1D0'
            write(stdout, "(A)") 'C11 = 232.4D0'
            write(stdout, "(A)") 'C22 = 119.36D0'
            write(stdout, "(A)") 'C44 = 45.95D0'
            write(stdout, "(A)") 'Bulk = 157.0D0'
            write(stdout, "(A)") '/'
            write(stdout, *)
            write(stdout, "(A)") "Don't forget that the last line must be empty!"
        end if
        read(fileid, nml=bound, iostat=iosuccess)
        if (iosuccess /= 0) then
            write(stdout, "(A)") "Error occuring during read ../parameter.nml, please check your format."
            write(stdout, "(A)") "The reference format of ../parameter.nml like:"
            write(stdout, "(A)") '&bound'
            write(stdout, "(A)") "re_bound = [0.0d0, 5.0d0]"       ! re参数的取值上下限
            write(stdout, "(A)") 'fe_bound = [0.0d0, 5.0d0]'       ! fe参数的取值上下限
            write(stdout, "(A)") 'rou_e_bound = [0.0d0, 50.0d0]'   ! rou_e参数的取值上下限
            write(stdout, "(A)") 'rou_s_bound = [0.0d0, 50.0d0]'   ! rou_s参数的取值上下限
            write(stdout, "(A)") 'alpha_bound = [0.0d0, 10.0d0]'   ! alpha参数的取值上下限
            write(stdout, "(A)") 'beta_bound = [0.0d0, 10.0]'      ! beta参数的取值上下限
            write(stdout, "(A)") 'A_bound =[0.0d0, 1.0d0]'         ! A参数的取值上下限
            write(stdout, "(A)") 'B_bound = [0.0d0, 1.0d0]'        ! B参数的取值上下限
            write(stdout, "(A)") 'cai_bound = [0.0d0, 1.0d0]'      ! cai参数的取值上下限
            write(stdout, "(A)") 'lamda_bound = [0.0d0, 1.0d0]'    ! lamda参数的取值上下限
            write(stdout, "(A)") 'Fn0_bound = [-5.0d0, 0.0d0]'     ! Fn0参数的取值上下限
            write(stdout, "(A)") 'Fn1_bound = [-2.0d0, 0.0d0]'     ! Fn1参数的取值上下限
            write(stdout, "(A)") 'Fn2_bound = [0.0d0, 2.0d0]'      ! Fn2参数的取值上下限
            write(stdout, "(A)") 'Fn3_bound = [-5.0d0, 0.0d0]'     ! Fn3参数的取值上下限
            write(stdout, "(A)") 'F0_bound = [-6.0d0, 0.0d0]'      ! F0参数的取值上下限
            write(stdout, "(A)") 'F1_bound = [0.0d0, 0.0d0]'       ! F1参数的取值上下限
            write(stdout, "(A)") 'F2_bound = [0.0d0, 3.5d0]'       ! F2参数的取值上下限
            write(stdout, "(A)") 'F3_bound = [-3.0d0, 1.0d0]'      ! F3参数的取值上下限
            write(stdout, "(A)") 'eta_bound = [-1.0d0, 2.0d0]'     ! eta参数的取值上下限
            write(stdout, "(A)") 'F_e_bound = [-5.0d0, 0.0d0]'     ! F_e参数的取值上下限
            write(stdout, "(A)") '/'
            write(stdout, *)
            write(stdout, "(A)") "Don't forget that the last line must be empty!"
        end if
        close(fileid)
    end subroutine read_target_value

    subroutine write_eam_file(lat)
        real(dp), intent(in) :: lat     !< 用于计算势函数的作用范围
        integer :: fileid, iosuccess
        logical :: file_exist
        real(dp) :: rc, rho_max
        real(dp) :: dr, drho
        real(dp) :: r, rou
        real(dp), allocatable, dimension(:) :: phi, rho, emb
        integer :: i
        character(len=8) :: date
        namelist/eamparameter/ re, fe, rou_e, rou_s, alpha, beta, A, B, cai, lamda, Fn0, Fn1, Fn2, Fn3, F0, F1, F2, F3, eta, F_e

        ! 如果是第一次, 那就从文件中读取数值
        if (.not. first_read_parameter) then
            inquire(file=filename, exist=file_exist)
            if (.not. file_exist) then
                write(stdout, "(A)") "File ../parameter.nml doesn't exist!"
                stop
            end if
            open(newunit=fileid, file=filename,action="read")
            read(fileid, nml=eamparameter, iostat=iosuccess)
            if (iosuccess /= 0) then
                write(stdout, "(A)") "Error occuring during read ../parameter.nml, please check your format."
                write(stdout, "(A)") "Don't forget that the last line must be empty!"
                stop
            end if
            close(fileid)
            first_read_parameter = .True.
        end if
        ! 计算rhol和rhoh
        rhoin = 0.85_dp * rou_e
        rhoout = 1.15_dp * rou_e
        ! 如果不是，就根据传入子程序的值写出对应eam文件
        rc = sqrt(10.0_dp) / 2 * lat            ! 计算势函数的截断半径
        rho_max = max(2.0d0*rou_e, 100.0_dp)    ! 计算势函数的最大电子密度
        dr = rc / (nr - 1)
        drho = rho_max / (nrho - 1)
        allocate(phi(nr), rho(nr), emb(nrho), source=0.0d0)
        ! 计算电子密度ρ(r)和两体势φ(r)
        do i = 1, nr
            r = (i-1) * dr
            if (r < 0.5d0) r = 0.5d0
            rho(i) = prof(r)
            phi(i) = pair(r)*r
        end do
        ! 计算嵌入能F(ρ)
        do i = 2, nrho
            rou = (i-1) * drho
            emb(i) = embed(rou)
        end do
        ! 将结果写入到文件中
        call date_and_time(DATE=date)
        open(newunit=fileid, file='../V.eam.alloy', action='write')
        ! 写入前三行的注释行内容
        write(fileid,*) 'DATE: ',date(1:4),'-',date(5:6),'-',date(7:8),   &
        ' UNITS: metal CONTRIBUTOR: Xiaowang Zhou xzhou@sandia.gov and ', &
        'Lucas Hale lucas.hale@nist.gov '
        write(fileid,*) 'CITATION: X. W. Zhou, R. A. Johnson, ', &
        'H. N. G. Wadley, Phys. Rev. B, 69, 144113(2004)'
        write(fileid,*) 'Generated by create.f and eam_fit'
        write(fileid,*) 1, 'V'
        write(fileid,*) nrho, drho, nr, dr, rc
        ! 写入V的相关信息
        write(fileid, '(3x, A)') "23        50.94000         3.03      bcc"
        ! 将嵌入能F(ρ)写入
        write(fileid, "(5es24.16)") (emb(i), i = 1, nrho)
        ! 将电子密度ρ(r)写入
        write(fileid, "(5es24.16)") (rho(i), i = 1, nr)
        ! 将两体势φ(r)写入
        write(fileid, "(5es24.16)") (phi(i), i = 1, nr)
        ! 关闭文件
        close(fileid)
        ! debug使用
        ! print *, rhoin
        ! print *, rhoout
    end subroutine write_eam_file

    ! 此函数用于计算电子密度
    pure function prof(r) result(f)
        real(dp), intent(in) :: r
        real(dp) :: f

        f = fe * exp(-beta*(r/re-1.0_dp))
        f = f / (1.0_dp+(r/re-lamda)**20)

    end function prof

    ! 此函数用于计算原子之间的两体势
    pure function pair(r) result(psi)
        real(dp), intent(in) ::r 
        real(dp) :: psi
        real(dp) :: psi1, psi2

        psi1 = ( A*exp(-alpha*(r/re-1)) ) / ( 1 + (r/re - cai)**20 )
        psi2 = ( B*exp(-beta*(r/re-1)) ) / ( 1 + (r/re - lamda)**20 )
        psi = psi1 - psi2
    end function pair 

    ! 此函数用于计算嵌入能F(ρ)
    pure function embed(rho) result(emb)
        real(dp), intent(in) :: rho
        real(dp) :: emb
        real(dp) :: dr 

        if (rho == 0.0d0) then
            emb = 0.0d0
        else if (rho < rhoin) then
            dr = rho / rhoin - 1
            emb = Fn0 + Fn1*dr + Fn2*dr**2 + Fn3*dr**3
        else if (rho < rhoout) then
            dr = rho / rou_e - 1
            emb = F0 + F1*dr + F2*dr**2 + F3*dr**3
        else 
            dr = rho/rou_s
            emb = F_e * (1.0_dp - eta*log(dr)) * dr**eta
        end if
    end function embed

end module file_IO