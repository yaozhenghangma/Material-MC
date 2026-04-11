# brainstorm: fix vesta kh bond vector

## Goal

修复 `src/spin_out.cpp` 中 KH (`Kitaev-Heisenberg`) 的 VESTA bond 向量导出，使其在周期性边界条件下输出正确的最短 bond 向量，并按用户确认的 `VECTR` 记录语义写出方向与起点原子索引。

## What I already know

* 目标函数是 `WriteVestaKhBondColor(Supercell&, std::string)`（`src/spin_out.cpp`）。
* 当前 `VECTR` 写法为两行：
  * 第一行：`vector_id delta_x delta_y delta_z 0`
  * 第二行：`vector_id tail_x tail_y tail_z 0`
* 当前实现的 `delta` 来自 wrapped 下标的笛卡尔坐标差 `target_cart - source_cart`，跨 PBC 时可能出现 ±1 个超胞晶格向量偏移。
* 当前去重通过 `visited_bonds`（无向 pair）保留首个遍历到的方向，方向可能受遍历顺序影响。
* 邻接在 `InitializeSupercell(...)` 中通过 modulo 包裹后物化为指针，未显式保存“跨边界平移量”。
* 用户已明确：PBC 应基于 `d = target_frac - source_frac` 后再做加减晶格常数（最小镜像）取最短向量。
* 用户已明确：`VECTR` 第二行应写“source 原子序号 + 四个 0”，而不是 tail 坐标。

## Assumptions (temporary)

* “source 原子序号”指 `STRUC` 段中的 1-based 原子序号（与输出原子顺序一致）。
* KH bond 颜色语义保持不变：x=红，y=绿，z=蓝（`VECTT` 逻辑保持兼容）。

## Open Questions

* （当前无阻塞问题）

## Requirements (evolving)

* 基于分数坐标计算 `d = target_frac - source_frac`。
* 对 `d` 做 PBC 最小镜像处理（按超胞周期 `n_x/n_y/n_z`），得到物理最短向量。
* `VECTR` 第一行中的 `d` 必须在“晶格方向归一化后的坐标系”下输出：
  * 取原 POSCAR 的晶格轴 `a/b/c`，分别归一化得到方向轴 `â/b̂/ĉ`。
  * 将最小镜像后的实际空间向量表示为该坐标系下的分量后写出。
* 将 source 原子序号写入 `VECTR` 第二行第一列，后四列写 0。
* 维持 KH bond 颜色输出（`VECTT`）不变。
* 保持输出可复现（相同输入下文件稳定）。
* `source_atom_index` 明确采用 `STRUC` 中当前输出顺序对应的 1-based 原子序号（i/j/k/l 四重循环递增编号）。

## Acceptance Criteria (evolving)

* [ ] 对跨 PBC 的 KH 邻接，导出的 bond 向量是最短镜像向量（无 ±1 晶格常数跳变）。
* [ ] `VECTR` 每个 bond 的第一行是更正后的 `d_x d_y d_z`，且分量定义在 `â/b̂/ĉ` 坐标系下。
* [ ] `VECTR` 每个 bond 的第二行首字段为 source 原子序号，其余为 0。
* [ ] `VECTT` 的 x/y/z 颜色映射仍为红/绿/蓝。
* [ ] 相同输入重复运行，`.vesta` 输出稳定一致。

## Definition of Done (team quality bar)

* Tests added/updated (unit/integration where appropriate)
* Lint / typecheck / CI green
* Docs/notes updated if behavior changes
* Rollout/rollback considered if risky

## Out of Scope (explicit)

* Hamiltonian 物理模型或 Monte Carlo 数值逻辑修改。
* 非 KH 模型的可视化行为变更。
* 除 `VECTR`/`VECTT` 相关记录外的大范围 VESTA 格式重构。

## Technical Notes

* 已检查文件：
  * `src/spin_out.cpp`（VESTA 输出实现）
  * `src/initialization.cpp`（neighbor 指针与 KH direction 模板生成）
  * `example/kh_minimal/VALIDATION_CHECKLIST.md`（KH 输出验证清单）
* 历史任务 `.trellis/tasks/archive/2026-04/04-11-fix-vesta-parallel-lines/prd.md` 中已记录同类问题背景，可复用其验证思路。

## Decision (ADR-lite)

**Context**: KH 的 VESTA 向量在跨 PBC 时出现错误长度/方向，同时当前第二行记录使用了 tail 坐标而不是用户期望的 source 原子索引。

**Decision**:
1. 采用分数坐标最小镜像策略：`d = target_frac - source_frac` 后按超胞周期做 wrap，取最短周期等价向量。
2. `VECTR` 输出采用用户确认语义：
   - 第一行 `vector_id d_x d_y d_z 0`，其中 `d_x/d_y/d_z` 是在归一化晶格方向基 `â/b̂/ĉ` 下的分量
   - 第二行 `source_atom_index 0 0 0 0`
3. `source_atom_index` 采用 `STRUC` 段当前输出顺序的 1-based 编号。

**Consequences**:
* 可消除跨边界 bond 的 ±1 晶格常数跳变伪影。
* 向量方向由 source→target 与最小镜像共同确定，行为可预测。
* 需要在导出阶段维护 site 到 `STRUC` 原子序号的确定性映射。

## Technical Approach

1. 在 `WriteVestaKhBondColor(...)` 中建立 `(i,j,k,l) -> atom_index`（1-based）映射（与 `STRUC` 写出顺序一致）。
2. 将 bond 向量计算切换到分数坐标差并应用最小镜像（按 `n_x/n_y/n_z` wrap）。
3. 由更正后的分数向量先转换到实际空间向量，再转换到 `â/b̂/ĉ` 坐标系分量并写入 `VECTR` 第一行。
4. `VECTR` 第二行写 source 原子序号与 0 填充。
5. 保持 `VECTT` 颜色逻辑不变。
