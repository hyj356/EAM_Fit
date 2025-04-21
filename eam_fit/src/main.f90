program main
    use iso_fortran_env, only: stdout => output_unit, dp => real64, sp => real32
    use file_IO, only: read_target_value, write_eam_file, generate_eam_parameter, summary, update_eam_parameter
    use optimize, only: cost_function, init_optimize
    implicit none
    integer, parameter :: max_iteration = 1000
    real(dp), parameter :: final_temperature = 1e-10_dp
    real(dp), parameter :: alpha = 0.99d0
    character(len=3) :: lat    !< 晶格类型
    real(dp) :: a0             !< 晶格常数
    real(dp) :: Ec             !< 内聚能
    real(dp) :: Efv            !< 空位形成能
    real(dp) :: C11            !< 弹性常数C11
    real(dp) :: C12            !< 弹性常数C22
    real(dp) :: C44            !< 弹性常数C44
    real(dp) :: Bulk           !< 体积模量
    real(dp) :: new_cost, old_cost, delta_cost
    real(dp) :: flag
    real(dp) :: initial_temperature = 100.d0, temperature
    integer  :: i
    
    ! 初始化随机数种子
    call random_seed()  
    ! 读取用户设置的拟合目标值, 以及对应拟合参数上下限
    call read_target_value(lat, a0, Ec, Efv, C11, C12, C44, Bulk)
    ! 写出eam/alloy格式文件
    call write_eam_file(a0)
    ! 记录初始的参数
    call update_eam_parameter(2)
    ! 初始化优化实例
    call init_optimize()
    ! 将old_cost设置为双精度浮点数的数值上限
    old_cost = huge(1.0d0)  
    ! 初始化退火温度
    temperature = initial_temperature
    ! 正式开始退火算法模拟迭代
    do i = 1, max_iteration
        ! 随机更新eam势函数参数
        call generate_eam_parameter()
        ! 将eam势函数参数写出
        call write_eam_file(a0)
        ! 调用lammps库计算当前势函数文件设置下与目标值的差距
        call cost_function(a0, Ec, Efv, C11, C12, C44, Bulk, new_cost)
        ! 计算与上一次的目标函数的值的差
        delta_cost = new_cost - old_cost
        ! 生成一个0-1之间均匀分布的随机数
        call random_number(flag)
        ! 利用metropelis准则, 判定是否接受新解, 
        if (delta_cost < 0.0d0 .or. flag < exp(-delta_cost / temperature)) then
            write(stdout, *) "Iteration: ", i, " accepted." , " cost: ", new_cost, &
                ' temperature: ', temperature
            old_cost = new_cost
            ! 如果接受新解, 记录这个新解
            call update_eam_parameter(2)
        else 
            write(stdout, *) "Iteration: ", i, " rejected.", " cost: ", new_cost, &
                ' temperature: ', temperature
            ! 如果不接收新解, 就回退到上一个解
            call update_eam_parameter(1)
        end if
        ! 降温
        temperature = temperature * alpha
        ! 如果温度达到了最低温度, 退出循环
        if (temperature < final_temperature ) exit
    end do
    ! 将eam势函数参数写出
    call write_eam_file(a0)
    ! 迭代结束, 打印总结
    call summary(old_cost)
    
end program main