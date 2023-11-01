
for(k in 2:15){
  for(f in 1:5){
    #submit job
    cmd <- paste('sbatch -N 1 --mem=10g -n 1 -t 08:00:00 --wrap="Rscript run_cv.R',
                 f,k,'"')
    system(cmd)
    
  }
}
