# SPECK
A Python adaptation of the R SPECK implementation in R. 

The following pseudocode was used as a guide, including the methods section described in the SPECK paper (). 


Algorithm 1 SPECK Algorithm

1: Xm,n−→−−−−−−NormalizationXm,n

Rank-k Selection and Singular Value Decomposition

2: ∑100i←1μiσivTi          ▹ Rank-100 SVD of Xm,n

3: rsdev←(σ1,σ2,…,σ100)2m−1−−√  ▹ PCs stdev.

4: rdiff←|rate of change(rsdev)|

5: rval←rle(rdiff).values      ▹ Run length encoding.

6: rlen←rle(rdiff).lengths

7: k←rdiff.index(min(rval≥0.01∧rlen≥2))+1      ▹ Rank-k.

8: Xm,n←μm,k×σk,k×vTk,n    ▹ Matrix reconstruction.

Clustered Thresholding

9: forv←1 to n do          ▹ Threshold each gene.

ckres←Ckmeans.1d.dp(x=Xm,n[v], k=c(1:4))

10:   num←ckres.cluster

11:   val←ckres.centers

12:   iflen(set(num))>1 then

13:    ind←val.index(min(val))+1

14:    Xm,n[v][num=ind]←0

15:   end if

16: end for
