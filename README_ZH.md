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
### 离线安装xmake
通过git下载源码到本地并上传到离线计算机。
```bash
git clone --recursive https://github.com/xmake-io/xmake.git
```
在源码所在文件夹进行编译安装。
```bash
make build
./scripts/get.sh __local__ __install_only__
source ~/.xmake/profile
```
xmake安装路径：`~/.local/bin`

通过git下载并更新xrepo官方仓库目录。

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
2. 参数文件input.txt

### 输出文件
1. 命令行输出
2. 物理学量统计值output.txt
3. 结构文件spin.txt