# Existence-of-non-convex-V-states
论文《不可压二维欧拉方程的 V-态存在性》的计算机辅助证明。

proof.cpp为引理2.6证明过程涉及的代码，如果存在$E[0](X_k)$的误差估计越界，程序会输出$X_k$及对应的$E[0](X_k)$误差范围，并最终输出“Proof is error!”，否则，程序输出“|E[x]| is bounded in 3e-8.”。

proof2.cpp为引理2.7证明过程涉及的代码，如果存在$E''[0](X_k)$的误差估计越界，程序会输出$X_k$及对应的$E''[0](X_k)$误差范围，并最终输出“Proof is error!”，否则，程序输出“|ddE[x]| is bounded in 50.”。

