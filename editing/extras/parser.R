# makeScripts: create N (.R + .sh) scripts for submitting to the cluster
makeScripts <- function(N = 1, origFile = NULL, folder = '', time = '0-01:00:00',
                        cpus = 1, msg = 1, mem = '2gb', jobName = NULL, breaks = NULL,
                        output = '~/Desktop', newFile = NULL, path = getwd(), ib = FALSE,
                        partition = 'sixhour', args = NULL, argtypes = TRUE,
                        modnets = FALSE, loop = TRUE){
  msgs <- c('FAIL', 'END,FAIL', 'START,END,FAIL')
  partition <- match.arg(tolower(partition), choices = c('crmda', 'sixhour'))
  if(!inherits(origFile, c('character', 'NULL'))){stop('Need to specify origFile')}
  if(!grepl('.R$', origFile)){origFile <- paste0(origFile, '.R')}
  if(is.null(jobName)){jobName <- ifelse(folder == '', 'X', gsub('^/|/$', '', folder))}
  if(is.null(newFile)){newFile <- jobName}
  list2env(lapply(list(folder = folder, output = output, path = path), function(z){
    paste0(gsub('/$', '', z), '/')}), envir = environment())
  if(folder == '/'){
    folder <- ''
  } else if(!gsub('/$', '', folder) %in% dir(output)){
    system(paste0('mkdir ', output, folder))
  }
  for(n in N){
    if(any(!is.null(breaks))){
      breaks <- breaks[order(breaks)]
      breaks <- breaks - c(0:(length(breaks)-1))
    }
    z <- as.character(parse(paste0(path, origFile)))
    z[1] <- paste0('n <- ', n)
    Z <- paste0(z[1], '\n', z[2])
    for(i in 3:length(z)){
      Z <- ifelse(i %in% breaks, paste0(Z, '\n\n', z[i]), 
                  paste0(Z, '\n', z[i]))
    }
    Rfile <- paste0(output, folder, newFile, n, '.R')
    cat(Z, file = Rfile)
    l <- vector('list', 17)
    l[1] <- '#!/bin/sh'
    l[2] <- '#'
    l[3] <- '#SBATCH --mail-user=t092s958@ku.edu'
    l[4] <- paste0('#SBATCH --job-name=', jobName, n)
    l[5] <- '#SBATCH --nodes=1'
    l[6] <- '#SBATCH --ntasks-per-node=1'
    l[7] <- '#SBATCH --constraint=ib'
    l[8] <- paste0('#SBATCH --time=', time)
    l[9] <- paste0('#SBATCH --mail-type=', msgs[ifelse(msg %in% 1:3, msg, 1)])
    l[10] <- paste0('#SBATCH --cpus-per-task=', cpus)
    l[11] <- paste0('#SBATCH --mem-per-cpu=', mem)
    l[12] <- paste0('#SBATCH --partition=', partition)
    l[13] <- '#SBATCH --error=\"%x.e%j\"'
    l[14] <- '#SBATCH --output=\"%x.o%j\"'
    l[15] <- 'cd $SLURM_SUBMIT_DIR'
    l[16] <- 'module load R'
    l[17] <- paste0('R --vanilla -f ~/', folder, newFile, n, '.R')
    if(!is.null(args)){
      if(is.character(args) & !is(args, 'list')){
        if(!grepl('=', args)){stop('Arg name and value must be separated by equals sign')}
        l[17] <- paste(l[17], '--args', gsub(' ', '', args))
      } else if(is(args, 'list')){
        if(is.null(names(args)) | any(names(args) == '')){stop('Must name all arguments')}
        if(!all(sapply(args, length) == 1) | any(sapply(args, mode) %in% c('list', 'function'))){
          stop('Each argument must contain a single element -- a string, number, or logical value')
        }
        args0 <- args
        if(isTRUE(argtypes)){
          names(args0) <- paste0(names(args0), toupper(sapply(
            sapply(args0, mode), substring, 1, 1)))
        }
        l[17] <- paste(l[17], '--args', paste(sapply(seq_along(args0), function(i){
          if(is.character(args0[[i]])){args0[[i]] <- paste0('"', args0[[i]], '"')}
          paste0(names(args0)[i], '=', args0[[i]])
        }), collapse = ' '))
      }
    }
    if(!ib){l <- l[-7]}
    ll <- length(l)
    Z <- paste0(paste0(l[1:(ll - 3)], collapse = '\n'), '\n\n', paste0(l[(ll - 2):ll], collapse = '\n\n'))
    shFile <- paste0(output, folder, newFile, n, '.sh')
    cat(Z, file = shFile)
  }
  if(isTRUE(modnets)){setupModnets(gsub('^/|/$', '', folder))}
  if(loop & length(N) > 1){
    loop <- c('#!/bin/sh\n', paste0('for i in ', paste(N, collapse = ' ')), 
              'do', paste0('  sbatch ', newFile, '$i.sh'), 'done')
    loop <- paste0(paste(loop, collapse = '\n'), '\n')
    cat(loop, file = paste0(output, folder, 'loop_submit.sh'))
  }
}

##### setupModnets: add 'modnets' folder to submission script(s) destination
setupModnets <- function(folders = '.', rm = FALSE){
  if(any(folders %in% c('.', ''))){folders <- '.'}
  fold <- "~/C:\\ Desktop/COMPS/METHODS/CODE/modnets/"
  mods <- paste0(c("functions", "ggm", "centrality", "sim", "mlGVAR", 
                   "simGVAR", "penalized", "power", "plots"), ".R")
  for(j in seq_along(folders)){
    if(!folders[j] %in% dir() & folders[j] != '.'){
      stopifnot(length(folders) == 1)
      folders[j] <- '.'
    }
    if(isTRUE(rm)){system(paste0('rm -rf ', folders[j], '/modnets'))}
    system(paste0('mkdir ', folders[j], '/modnets'))
    for(k in seq_along(mods)){system(paste0("cp ", fold, mods[k], " ./", folders[j], "/modnets/"))}
    f0 <- readLines(paste0(folders[j], "/modnets/functions.R"))
    f1 <- which(f0 == "setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')")
    f2 <- which(f0 == "files <- paste0('./', c('ggm', 'centrality', 'sim', 'mlGVAR', 'simGVAR', 'penalized', 'power', 'plots'),'.R')")
    f0[f1] <- paste0("#", f0[f1])
    f0[f2] <- gsub("'./'", "'modnets/'", f0[f2])
    f0 <- paste0(paste0(f0, collapse = "\n"), "\n")
    cat(f0, file = paste0(folders[j], "/modnets/functions.R"))
  }
}
