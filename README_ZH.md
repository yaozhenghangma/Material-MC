# Material-MC

[![build](https://github.com/yaozhenghangma/Material-MC/actions/workflows/main.yml/badge.svg)](https://github.com/yaozhenghangma/Material-MC)
[![release](https://img.shields.io/github/release/yaozhenghangma/Material-MC)](https://github.com/yaozhenghangma/Material-MC/releases)
[![Common Changelog](https://common-changelog.org/badge.svg)](https://github.com/yaozhenghangma/Material-MC/blob/main/CHANGELOG.md)

晶体通用蒙特卡洛模拟

[**中文版**](./README_ZH.md) ｜ [**English**](./README.md)

## 编译
### 依赖项
- xmake
- MPI
- Boost.mpi
- Boost.serialization
### 

## 使用方法
### 命令行输入
帮助文档
```bash
./MMC -h
```
64进程并行运行程序
```bash
mpirun -np 64 ./MMC
```

### 输入文件
1. 结构文件POSCAR
2. 参数文件input.toml

### 输出文件
1. 日志文件log.txt
2. 物理学量统计值output.txt
3. 结构文件spin.xsf