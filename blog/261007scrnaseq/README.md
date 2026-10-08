# AWS 上比较 Cell Ranger 与 Simpleaf：CloudFormation 模板

模板创建一台临时 Linux 服务器，自动安装软件、先跑小测试，再按你选择的模式运行完整实验。无需在 Mac 安装 Docker，也无需 SSH 密钥。

**2026-10-05 已在 AWS us-east-1 实测：tests-only 与 full 均成功完成，结果在 `scrna-results/`。实测记录见第 8 节。**

## 1. 上传哪个文件

只上传 `scrnaseq-comparison.yaml`。它已经包含安装脚本和运行脚本，不需要另找脚本文件。

在 AWS 控制台选择 **US East (N. Virginia)，us-east-1**，打开 CloudFormation → Create stack → With new resources → Upload a template file。

1. 上传 YAML，填写 Stack name，例如 `scrna-test`。
2. 首次保持 `RunMode = tests-only`，磁盘 `DiskGiB = 250`，时间限制可选 `MaxRuntimeHours = 4`。`SaveAlignments` 和 `SuccessAction` 保持默认（见第 2、6 节）。
3. 其他选项保持默认。最后勾选允许 CloudFormation 创建 IAM 资源，然后创建。
4. 在 Outputs 里找到 `ResultsBucketName`，到对应 S3 桶查看 `runs/<时间>/metadata/status.json`。
5. 测试状态为 `COMPLETED` 后，删除测试栈。默认 `SuccessAction = terminate` 时，EC2 在上传核对成功后已自行终止；删除栈释放剩余的网络和 IAM 资源。测试结果的 S3 桶会保留。
6. 用同一 YAML **创建一个新栈**，例如 `scrna-full`，设置 `RunMode = full`、`MaxRuntimeHours = 24`，其余不变。

`full` 也会先跑小测试，成功后才下载完整数据。不要通过更新测试栈来切换模式：修改 UserData 不保证自动重新执行。

**CloudFormation 的 `CREATE_COMPLETE` 表示服务器环境初始化成功，不表示分析完成。** 分析成功以 S3 的 `status.json` 中 `state = COMPLETED` 为准。

## 2. 自动做什么

顺序为：创建网络和服务器 → 安装 Java、Docker、Nextflow → Cell Ranger 小测试 → Simpleaf 小测试 → 下载完整数据与参考 → 核对文件并生成样本表 → Cell Ranger → Simpleaf＋QCatch → 上传结果并核对 → 成功则终止服务器，失败则停机。

两条路线顺序执行，任务也采用保守的串行调度，不是性能竞赛。Nextflow 的时间报告包含索引构建和其他辅助步骤；不能直接把整个流程耗时当成单个工具速度。

| 模板设置 | 默认值与用途 |
| --- | --- |
| 服务器 | `r6i.4xlarge`，x86_64，16 vCPU、128 GiB RAM，按需实例 |
| 操作系统 | Canonical Ubuntu 22.04 LTS；AMI 由区域公共 SSM 参数解析 |
| 磁盘 | 加密 250 GiB gp3；`/data/scrna` 与系统共用此磁盘。估计峰值 110–150 GB；解压后删除两个压缩包；每分钟记录已用空间到 `metadata/disk-used-gib.log` |
| 分析资源 | 最多 12 CPU、112 GB；主要分析步骤 96 GB，给系统和 Nextflow 留空间 |
| 网络 | 独立 VPC、公有子网、公网 IPv4、互联网网关；无 NAT 网关，无入站端口 |
| 登录 | SSM Session Manager；不需要开放 SSH 或保存 AWS 密钥到服务器 |
| 结果 | 非公开、加密 S3 桶；删除栈时保留 |
| 自动收尾 | 上传并核对成功后**终止**实例（磁盘随之删除）；失败、上传未核对通过或 `SuccessAction = stop` 时**停机**保留磁盘；另有最长运行时间定时器（到时停机） |

实例启动需要该区域有足够的 EC2 按需 vCPU 配额和 `r6i.4xlarge` 容量。模板不会调整你的账户配额。

## 3. 完整实验下载什么

**数据：10x PBMC 1k v3，健康人外周血单核细胞，10x 3′ v3。** 使用整个公开 FASTQ 包，不抽样。

- FASTQ：<https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar>
- 数据介绍：<https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-3-standard-3-0-0>
- 共同参考：<https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz>
- 参考说明：<https://www.10xgenomics.com/support/software/cell-ranger/downloads>

Cell Ranger 直接使用这个参考包中的预建索引；Simpleaf 使用包中**相同的 FASTA/GTF**创建自己的索引。两种工具的索引格式和方法不同。

参考压缩包校验官方 MD5：`a7b5b7ceefe10e435719edc1a8b8b2fa`。下载后记录两个压缩包的 SHA256 和参考 FASTA/GTF 的 SHA256；FASTQ 逐个检查 gzip 完整性。记录 SHA256 是为了追溯文件，并不代表 FASTQ 已与官方校验值比对。

样本表自动配对所有 `_R1_` / `_R2_` FASTQ，缺少配对会报错。所有 lane 使用同一个样本名 `PBMC1K`，不将 I1 文件当表达序列。

没有填 `expected_cells`，也没有强制保留 1000 个细胞。公开数据的历史细胞数不是本实验必须达到的数值。这里使用更新的 GRCh38-2024-A 参考，因此也不是对原始 10x 发布结果的精确复刻。

## 4. 固定版本与比较边界

- nf-core/scrnaseq：`4.2.0`，检验源码 commit `3fc17b4f971a89e47c88337de71d0e777ffad8cc`。
- Nextflow：`25.10.4`。
- 使用该流程定义的容器；保存运行软件版本和本地 Docker 镜像摘要。OS 安全更新、容器注册表和插件下载仍依赖外部服务，不等于全部依赖已离线冻结。
- 分析运行时不加载社区集群配置；本机用 `local` executor＋Docker。单独保存的 `resolved.config` 是配置预览，不包含运行时参数，不能当成完整分析参数文件。
- 两条路线都关闭 CellBender。
- 完整运行默认 `SaveAlignments = true`：保存 Cell Ranger BAM（`--create-bam true`）和 Simpleaf `af_map`，用于在 read 层面追查有争议的 barcode，S3 约多 5–10 GB。设为 `false` 则不保存。小测试始终不保存。
- Simpleaf 使用 `cr-like` UMI resolution，开启 QCatch，关闭双细胞移除。

**本实验比较 Cell Ranger 的原生筛选结果与 Simpleaf＋QCatch 的结果。** QCatch 的筛选包含质控处理，不能将最终细胞数差异全部归因于比对器。Simpleaf 的计数结果也经过自身 barcode permit-list 处理，并不自动等同于包含所有液滴的 Cell Ranger raw matrix。

保留原始发布矩阵、筛选矩阵和各阶段报告。后续在 Mac 分析时，再核对 barcode、矩阵范围及 QC 标准；统一细胞识别方法需要另做，模板没有声称完成这个下游实验。

小测试使用流程提供的小鼠 chr19 数据，只验证执行链路，不参与 PBMC 比较；测试 profile 默认跳过 QCatch，完整模式才运行 QCatch。Cell Ranger 小测试生成的 `cellranger/mkref`（小鼠 chr19 索引，约 1 GB）不上传 S3。

## 5. 去哪里看进度与结果

S3 桶的结构：

```text
runs/<UTC时间>/
  results/          # 两条路线的发布矩阵、报告、软件版本及小测试输出
  metadata/         # 参数、样本表、运行命令、校验值、Nextflow 报告和状态
  task-logs/        # work 里的 .command.* 执行记录
bootstrap/          # 初始化失败时尽力上传的安装日志
```

`status.json` 的状态：

| state | 含义 |
| --- | --- |
| `RUNNING` | 正在运行；`phase` 指出最近的阶段 |
| `COMPLETED` | 两个小测试成功；若为 full，完整两条路线也成功，最终上传成功并经 dry-run 核对 |
| `FAILED` | 安装后的运行脚本失败，上传了能取得的日志和已有结果 |
| `UPLOAD_FAILED` | 最终上传失败或核对未通过；S3 可能不完整，实例停机保留磁盘，应先检查 EC2 磁盘 |

`status.json` 中的 `instance_action` 记录收尾方式：`terminate` 或 `stop`。

状态在阶段之间更新，不是实时心跳。若触发最长时间限制，可能还显示 `RUNNING`：服务器会直接关机，未保证最后一次上传。

运行期间，在 EC2 控制台选实例 → Connect → Session Manager，可以输入：

```bash
tail -n 60 /data/scrna/metadata/runner.log
sudo journalctl -u scrna-experiment.service -n 60 --no-pager
cat /data/scrna/metadata/status.json
df -h /data
```

Session Manager 页面显示尚未在线时，等待初始化；初始化故障先查看 EC2 的系统日志及 CloudFormation Events。

若本地已经配置 AWS CLI，可把 S3 下载到 Mac。将 `YOUR_BUCKET` 替换为 Outputs 中的名字：

```bash
# 只下载分析结果：跳过 BAM、Simpleaf af_map 和小测试输出（本地空间有限，它们留在 S3）
aws s3 sync s3://YOUR_BUCKET/runs/ ./scrna-results/ --region us-east-1 \
  --exclude '*.bam' --exclude '*.bam.bai' --exclude '*/af_map/*' \
  --exclude '*/results/test-cellranger/*' --exclude '*/results/test-simpleaf/*'
```

## 6. 失败、重启与费用

成功（上传后 dry-run 核对 `results/` 与 `task-logs/` 无遗漏）时 EC2 **自行终止**，磁盘随之删除，之后只剩 S3 存储费用。终止权限限定为本栈创建的实例。终止请求失败时退回为停机。

失败、上传未核对通过，或 `SuccessAction = stop` 时 EC2 **停止而非删除**：计算计费结束，但 EBS（250 GiB 约每天 0.65 美元）和 S3 仍收费。无论哪种情况都不会删除 S3 中的分析结果。

实例自行终止后，栈里的 EC2 资源已不存在：之后只做删除栈，不要更新栈。

最长时间从初始化脚本启动算起，包含安装、下载和测试；每次重启服务器定时器会重新开始。它是时间上限，不是金额上限。若任务较慢，24 小时可能不够，先查看日志再决定是否调整。

**运行中下载的 FASTQ、索引和完整 work 目录只留在 EC2 磁盘，不上传 S3。** 因此 S3 结果包不是完整的 `-resume` 备份；删除栈会删除这些本地文件。

初始化成功后，运行服务没有设置开机自动重跑，重启 EC2 可以先查看问题。不要直接重新运行 `run.sh`：它会重新克隆/下载，不是恢复脚本。需要恢复时保留原工作目录，在 `/data/scrna` 下复制 `metadata/commands.log` 中失败的命令，增加 `-resume`，确认原来的 `-work-dir` 和输入未改变；修复后手动上传结果并关机。模板未自动实现这一步。

```bash
cd /data/scrna
export NXF_HOME=/data/scrna/nextflow-home
export NXF_VER=25.10.4
export NXF_TEMP=/data/scrna/tmp TMPDIR=/data/scrna/tmp
export NXF_OPTS='-Xms512m -Xmx4g'
# 检查 metadata/commands.log；在失败命令末尾添加 -resume 后执行。
# 恢复成功后上传已有结果，最后主动关机。
sudo shutdown -h now
```

下载并确认需要的矩阵和日志后，删除 CloudFormation 栈，释放 EC2、EBS 和网络资源。S3 桶按设计保留；不需要时单独清空并删除它。不要在失败日志或重要本地结果尚未取回时删除栈。

Cell Ranger 及其他软件仍受各自许可条款约束；使用流程容器不改变软件许可。

## 7. 已验证与尚未验证

本地检查 CloudFormation 资源定义、UserData 大小、Bash 语法、内嵌 Python 语法、参数名、FASTQ 配对的正常/异常案例及上传失败处理。详见 `VALIDATION.txt`。

2026-10-05 修改后（磁盘 250 GiB、成功后自行终止、跳过测试 mkref 上传、保存比对文件）另做的本地检查：cfn-lint 1.57.1 无报错；UserData 渲染后 12.1 KB（EC2 上限按编码前 16 KB 计）；Bash 语法通过；用模拟的 `aws`/`curl` 测试了五种收尾情况：成功、分析失败、上传失败、核对发现遗漏、终止请求失败。这些检查都不代替在 AWS 上实测。

2026-10-05 已在 AWS 上实测通过，见第 8 节。在其他账户或区域使用时，仍需确认配额和所选可用区是否提供 `r6i.4xlarge`；首次使用先运行 tests-only。

参考文档：

- <https://nf-co.re/scrnaseq/4.2.0/docs/usage/>
- <https://nf-co.re/scrnaseq/4.2.0/docs/output/>
- <https://github.com/nf-core/scrnaseq/tree/4.2.0>
- <https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-attribute-creationpolicy.html>
- <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html>
- <https://docs.aws.amazon.com/systems-manager/latest/userguide/session-manager.html>

**English recap:** Upload the YAML to CloudFormation, run tests-only first, then create a new full stack. The experiment compares Cell Ranger with Simpleaf plus QCatch on the same PBMC FASTQs and reference annotation. Results go to retained S3; the instance stops automatically. Local validation passed; no AWS deployment or real-data run has been performed.

## 8. 实测记录（2026-10-05，us-east-1a）

成功的完整运行：栈 `scrna-full2`，run `20261005T210235Z`，从 21:02 运行到 22:26 UTC（约 84 分钟，含安装、两个小测试、下载）。上传核对通过后实例自行终止。

| 步骤 | 墙钟 | CPU 利用 |
| --- | --- | --- |
| 两个小测试 | 约 10 分钟 | |
| 下载 FASTQ（5.3 GB）+ 参考（10.6 GB）并校验 | 约 7 分钟 | |
| `CELLRANGER_COUNT`（预建索引，含 BAM） | 22 分 54 秒 | 641% |
| `SIMPLEAF_INDEX`（GRCh38 splici） | 16 分 18 秒 | 630% |
| `SIMPLEAF_QUANT` | 9 分 46 秒 | 513% |
| `QCATCH` | 3 分 29 秒 | 113% |

- 磁盘峰值 66 GiB（`metadata/disk-used-gib.log`），250 GiB 绰绰有余。
- S3 结果共 8.3 GiB / 867 个对象，其中大部分是 Cell Ranger BAM 和 Simpleaf `af_map`。
- 版本：`quay.io/nf-core/cellranger:10.0.0`，`quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0`，QCatch 来自 `community.wave.seqera.io/library/pip_qcatch`；镜像摘要见 `metadata/docker-images.txt`。
- 初步数字：Cell Ranger 估计 1,221 个细胞（Median UMI 10,029，内含子 reads 31.1%）；Simpleaf＋QCatch 保留 1,228 个细胞（Median UMI 10,531）。数量接近不代表 barcode 相同，重叠情况待在 Mac 上比较。

实测中修复的三个模板问题（之前的失败运行都在这些步骤停止，按设计停机保留了日志）：

1. `--custom_config_base ''`：Nextflow 把空字符串解析为 `true`，去找 `true/nfcore_custom.config`。改为指向本地空配置目录 `/data/scrna/empty-configs`。
2. 参数校验时访问 iGenomes 默认地址 `s3://ngi-igenomes/` 返回 403（使用实例角色签名）。改为 `--igenomes_ignore true --igenomes_base /data/scrna/empty-configs`。
3. `nextflow config ... -c file`：Nextflow 25.10 要求 `-c` 写在子命令之前。已改正，且该预览失败不再中断分析。

另外在子网上固定了可用区（默认 `us-east-1a`）：这个账户的 us-east-1e 不提供 `r6i.4xlarge`。

费用估算（按需价格，以 Cost Explorer 1–2 天后的账单为准）：五次实例共运行约 2.2 小时，r6i.4xlarge 约 1.01 美元/小时，计算约 2.2 美元；EBS、公网 IPv4 合计不到 0.1 美元；S3 存储 8.3 GiB 约每月 0.2 美元；下载到 Mac 的 8.3 GiB 若在每月 100 GB 免费流出额度内则不收费。

