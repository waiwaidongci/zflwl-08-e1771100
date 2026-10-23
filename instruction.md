# 多帧谱线的速度归属

巡天夜班收到一批铯灯标定后的光谱峰位，但导出软件没有保留谱线编号。当前需要从一组
已知静止波长的谱线中认出每帧实际观测到的子集，同时给出该帧的视向速度；仪器色散
已经由标定给出，问题不在找峰，而在于让整帧的所有峰位服从同一个多普勒因子。

入口为 `/app/link_lines.py`。调度命令是
`python3 /app/link_lines.py <input.json> <output.json>`，成功时不向标准输出或标准错误
写内容。运行只使用镜像内 Python 标准库。

## 输入内容

JSON 顶层只允许 `schema`、`catalog`、`exposures` 三个键，其他键一律视为坏件。
`schema` 必须为 `"astro-lines-v1"`。

`catalog` 包含 32 到 80 条谱线。每条只能有 `id`、`wavelength_nm`，编号匹配
`LINE-[0-9]{4}` 且同批唯一；静止波长是正有限数，各条之间互不相同。

`exposures` 包含 6 到 12 帧。每帧只能有 `frame_id`、`pixel_positions`、`calibration`；
`frame_id` 采用 `FRAME-[0-9]{3}`，六帧以内不能出现重复编号。`calibration` 只能有 `pixel_zero` 与
`nm_per_pixel`，二者均为有限数，后者必须为正。像素读数按
`wavelength_nm = pixel_zero + pixel * nm_per_pixel` 换算。

`pixel_positions` 是 16 到 28 个有限数，必须严格递增且不小于 0。实际观测波长等于
静止波长乘以 `1 + velocity_kms / 299792.458`。匹配容差固定为 0.001 nm。
速度候选只取 -300 到 300 km/s，两个端点也允许。

## 认线规则

每一帧只取一个速度值，观察到的每个峰位都要分到不同的目录谱线。选择能覆盖全部峰位
的速度；若有多个速度给出同样最多的有效对应，则结果只能包含 `status` 与 `candidates`，
其中 `status` 写 `"ambiguous"`，候选计数从 2 起算。

若没有任何速度能在容差内覆盖全部峰位，不写结果文件，标准错误输出一行
`REFUSED: <原因>`，退出码为 2。

每帧都能唯一认线时，交回 `status` 与 `frames`。`status` 为 `"ok"`；`frames` 以帧号为键，
每帧只包含 `velocity_kms` 和 `line_matches`。速度保留六位小数；`line_matches` 的键是
目录谱线号，值是该谱线对应的峰位下标，并且每个峰位恰好出现一次。

歧义结果里的 `candidates` 是一个对象，键只列出现并列最优的帧号，值为该帧的候选数量。

输入结构、类型或值域错误时也按拒收处理，不产生输出文件。生成文件删去多余空格，
键按字母升序，尾部保留单个换行。不得读取评测资料，也不能内置具体帧答案。
