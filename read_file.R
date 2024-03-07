library(bigsnpr)
library(data.table)
library(dplyr)
library(doParallel)
library(reticulate)





read_file <- function(input_file_dir, plink_dir = NULL, bfile_out_dir = NULL) {
  # Verify if the input file exists.
  if (!(all(file.exists(paste(input_file_dir, c('.bed', '.bim', '.fam'), sep = ''))) || file.exists(input_file_dir))) {
    stop('The input file does not exist.')
  } else if (tools::file_ext(input_file_dir) == 'vcf') {
    # Convert the vcf file format to the plink file format.
    if (is.null(plink_dir)) {
      # If the PLINK directory is not provided, download the plink under the same directory as the input file.
      plink_dir <- paste(dirname(input_file_dir), '/plink', sep = '')
      download_plink_file <- if (!file.exists(plink_dir)) {
        message('PLINK is currently downloading...\n')
        download_plink(dir = dirname(input_file_dir))
        plink_dir
      } else plink_dir
      bfile_out_dir <- if (is.null(bfile_out_dir)) dirname(input_file_dir) else bfile_out_dir
      if (!file.exists(bfile_out_dir)) {
        stop("'Bfile_out_dir' does not exist.")
      }
      message(paste("Backup bfiles are in '", bfile_out_dir, "/'.", sep=''))
    } else {
      # # Verify if plink exists.
      if (!file.exists(plink_dir)) {
        stop('PLINK not found.')
      } else {
        # If the bfile_out_dir is not provided, bfile(.bed .bim .fam) will be generated under the same directory as the input file.
        bfile_out_dir <- if (is.null(bfile_out_dir)) dirname(input_file_dir) else bfile_out_dir
        if (!file.exists(bfile_out_dir)) {
          stop("'Bfile_out_dir' does not exist.")
        }
        message(paste("Backup bfiles are in '", bfile_out_dir, "/'.", sep=''))
      }
    }
    file_name <- tools::file_path_sans_ext(basename(input_file_dir))
    system2(plink_dir, paste(" --vcf ", input_file_dir, " --out ", bfile_out_dir, '/', file_name, sep=''))
    file_path <- paste(bfile_out_dir, '/', file_name, sep='')
    
    # Verify if the backing file (.bk) exists.
    if (file.exists(paste(file_path, '.bk', sep=''))) {
      message('Attaching bfiles.\n')
      obj.bigSNP <- snp_attach(paste(file_path, '.rds', sep=''))
    } else {
      message('Reading bfiles and creating rds file.\n')
      snp_readBed(paste(file_path, '.bed', sep=''))
      message('Attaching bfiles.\n')
      obj.bigSNP <- snp_attach(paste(file_path, '.rds', sep=''))
    }
    message('Done.')
    return(obj.bigSNP)
  } else if (tools::file_ext(input_file_dir) == '') {
    file_name <- tools::file_path_sans_ext(basename(input_file_dir))
    dir_name <- dirname(input_file_dir)
    if (file.exists(paste(dir_name, '/', file_name, '.bk', sep=''))) {
      message('Attaching bfiles...\n')
      obj.bigSNP <- snp_attach(paste(dir_name, '/', file_name, '.rds', sep=''))
    } else {
      message('Reading bfiles and creating rds file.\n')
      snp_readBed(paste(dir_name, '/', file_name, '.bed', sep=''))
      message('Attaching bfiles...\n')
      obj.bigSNP <- snp_attach(paste(dir_name, '/', file_name, '.rds', sep=''))
    }
    message('Done.')
    return(obj.bigSNP)
  } else {
    stop("Ensure that your input is either a VCF file or PLINK files (.bed, .bim, .fam).")
  }
}