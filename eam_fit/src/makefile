# 编译器
FC = gfortran
# 编译选项
FFLAGS = -O3 -Wno-unused-value -fcheck=all
FFLIB  = -L/home/hyj356/soft/lammps-29Aug2024/src -llammps

# 获取当前目录下所有的 .f90 文件
# 生成对应的 .o 文件列表，存放在 obj 目录下
OBJDIR = obj
OBJS = lammps.o mod_FileIO.o mod_optimize.o main.o
# 可执行文件名称
EXEC = eamfit
# 默认目标
all: $(EXEC)

# 链接目标文件生成可执行文件
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS)  -o $@ $^ ${FFLIB}

# 编译 .f90 文件生成 .o 文件到 obj 目录
%.o: %.f90
	$(FC) $(FFLAGS) ${FFLIB} -c $< -o $@

# 清理生成的文件
clean:
	rm -f $(OBJS) $(EXEC)
	rm *.mod