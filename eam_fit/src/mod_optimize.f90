module optimize
    use iso_fortran_env, only: stdout => output_unit, dp => real64, sp => real32
    use iso_c_binding, only: c_double
    use LIBLAMMPS
    implicit none
    
    private
    public :: cost_function, init_optimize
    real(dp), parameter, private :: stress_factor = 0.1_dp  ! 弹性常数的误差乘以的缩放系数
    type(lammps), private :: lmp !< 一个包含lammps的实例
    CHARACTER(LEN=12), PARAMETER :: args(3) = [ CHARACTER(LEN=12) :: 'liblammps', '-screen', 'none' ]

contains

    ! 初始化lammps实例, 调用命令行指令使得lammps不会输出任何信息到屏幕上
    subroutine init_optimize()
        lmp = lammps(args)   
    end subroutine

    subroutine cost_function(a0, Ec, Efv, C11, C12, C44, Bulk, cost)
        real(dp), intent(in) :: a0             !< 期望的晶格常数
        real(dp), intent(in) :: Ec             !< 期望的内聚能
        real(dp), intent(in) :: Efv            !< 期望的空位形成能
        real(dp), intent(in) :: C11            !< 期望的弹性常数C11
        real(dp), intent(in) :: C12            !< 期望的弹性常数C22
        real(dp), intent(in) :: C44            !< 期望的弹性常数C44
        real(dp), intent(in) :: Bulk           !< 期望的体积模量
        real(dp), intent(inout) :: cost        !< 计算当前目标函数的函数值
        real(dp) :: a0_T                 !< 根据当前势函数参数设置计算出来的实际晶格常数
        real(dp) :: Ec_T                 !< 根据当前势函数参数设置计算出来的实际内聚能
        real(dp) :: Ev_T
        real(dp) :: C11_T, C12_T, C44_T, Bulk_T

        ! 执行晶格常数和内聚能测试
        call lmp%file('../lattice_cohesive/test.lmp')
        a0_T = lmp%extract_variable("lat")
        Ec_T = lmp%extract_variable("coh")
        ! 执行缺陷形成能测试
        call lmp%file("../vacant_formation/test.lmp")
        Ev_T = lmp%extract_variable("Ev")
        ! 执行弹性常数测试
        call lmp%file("../elastic_constant/in.elastic")
        C11_T = lmp%extract_variable("C11all")
        C12_T = lmp%extract_variable("C12all")
        C44_T = lmp%extract_variable("C44all")
        Bulk_T = lmp%extract_variable('bulkmodulus')
        ! 评估当前势函数输入下, 与目标值的差距
        cost = abs(a0 - a0_T) + abs(Ec - Ec_T) + abs(Efv - Ev_T) + &
               stress_factor * abs(C11 - C11_T) +                  &
               stress_factor * abs(C12 - C12_T) +                  &
               stress_factor * abs(C44 - C44_T) +                  &
               stress_factor * abs(Bulk - Bulk_T)
        ! 将结果输出
        ! write(stdout, *) a0_T, Ec_T, Ev_T
        ! write(stdout, *) C11_T, C12_T, C44_T  
    end subroutine cost_function

end module optimize