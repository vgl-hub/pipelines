insert in your ~/.bashrc or ~/.bash_profile

```
alias sq='squeue -p vgl,vgl_bigmem,vglbfx -o "%.15i %.15u %.15P %.15R %.15T %.15M %.j" | tr -s \\t'
alias sb='sbatch -p vgl -n 1 -c 32'
alias snodes='sinfo -N -p vgl,vglbfx,vgl_bigmem'
alias snodesgpu='sinfo -N -p hpc_a10,hpc_a100,hpc_k80,hpc_v100'
alias "l=ls -lahF"
## set prompt to have current host and wd
export PS1="[\u@\h \w]\$ "
export SACCT_FORMAT="JobID%20,JobName,User,Partition,NodeList,Elapsed,State,ExitCode,MaxRSS,AllocTRES%32"
```

examples of usage:
```
(VGP) [labueg@login04 ~]$ sq
          JOBID            USER       PARTITION NODELIST(REASON)           STATE            TIME NAME
       45705592      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g31460_bwa_mem2_labueg_rockefeller_edu
       45705593      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g31461_bwa_mem2_labueg_rockefeller_edu
       45729448      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g32833_bwa_mem2_labueg_rockefeller_edu
       45731656      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g32868_bwa_mem2_labueg_rockefeller_edu
       45731662      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g32867_bwa_mem2_labueg_rockefeller_edu
       45739793      vgl_galaxy      vgl_bigmem     (Resources)         PENDING            0:00 g32869_bwa_mem2_labueg_rockefeller_edu
       45677343      vgl_galaxy             vgl         node147         RUNNING      7-22:09:41 g32273_hifiasm_labueg_rockefeller_edu
       45731316   smrtanalysis2             vgl         node144         RUNNING      4-16:00:59 cromwell_f4f0d490_lima
       45752100          labueg             vgl         node143         RUNNING      2-13:29:33 aEleCoq3
       45756089       gformenti             vgl         node145         RUNNING      1-19:58:45 novoplasty.sh
       45760307       smarcus01             vgl         node144         RUNNING      1-07:31:01 sys/dashboard/sys/rstudio

(VGP) [labueg@login04 ~]$ snodes
NODELIST   NODES  PARTITION STATE
node141        1     vglbfx idle
node142        1     vglbfx idle
node143        1        vgl alloc
node144        1        vgl alloc
node145        1        vgl alloc
node146        1        vgl idle
node147        1        vgl alloc
node148        1        vgl idle
node149        1        vgl drain
node150        1        vgl idle
node151        1        vgl idle
node152        1        vgl idle
node153        1        vgl idle
node154        1        vgl idle
node155        1        vgl idle
node156        1        vgl idle
node157        1        vgl idle
node158        1        vgl idle
node159        1        vgl idle
node160        1        vgl idle
node161        1        vgl idle
node162        1        vgl idle
node163        1        vgl idle
node164        1        vgl idle
node165        1        vgl idle
node166        1        vgl idle
node167        1        vgl idle
node168        1        vgl idle
node169        1        vgl idle
node170        1        vgl idle
node171        1 vgl_bigmem drain

(VGP) [labueg@login04 ~]$ sb --wrap="hostname"
Submitted batch job 45767980
(VGP) [labueg@login04 ~]$ cat slurm-45767980.out
node146.hpc.rockefeller.internal
```