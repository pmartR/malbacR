context('Run optimal_span which is a get_params (which is a subsection of bc_qcrlsc)')

test_that('optimal_span returns the optimal span for QCRLSC Normalization', {
  
  # Load the reduced metabolite data frames ------------------------------------
  load(system.file('testdata',
                   'opt_span.rda',
                   package = 'malbacR'))

  # run through the checks -----------------------------------------------------
  # data must be a data frame
  expect_error(optimal_span(opt_span,lam = c(1,2)),
               "qc_data must be a data.frame")
  # make the data data.frame
  data2run <- opt_span[[1]]
  
  # data must be two columns (here we have 3)
  expect_error(optimal_span(data2run,lam = c(1,2)),
               "qc_data must have 2 columns")
  # remove the batch column
  data2run <- data2run[,-3]
  
  # check that lam is a numeric vector
  expect_error(optimal_span(data2run,lam = c("dog","cat")),
               "lam must be a numeric vector")
  expect_error(optimal_span(data2run,lam = list(1,2)),
               "lam must be a numeric vector")
  
  # Run optimal spans ----------------------------------------------------------
  opt_output <- optimal_span(data2run,lam = c(1,2))
  
  # Dimension check time -------------------------------------------------------
  # we should have  a list of length 4
  expect_equal(names(opt_output),c("Lambda","Alpha","MSE","all_res"))
  
  # lambda should be a single number
  expect_equal(length(opt_output$Lambda),1)
  # alpha should be a single number
  expect_equal(length(opt_output$Alpha),1)
  # mse should be a single number
  expect_equal(length(opt_output$MSE),1)
  
  # all_res should be a data.frame with 3 columns
  expect_equal(colnames(opt_output$all_res), c("Lambda","Alpha","MSE"))
  # nrow is based on sequence of possible alpha values that are considered
  
  # find the potential alpha values
  lam = c(1,2)
  alp_seq1 = seq(from = (lam[1] +1), to = (nrow(data2run)-1), by = 1)/(nrow(data2run)-1)
  alp_seq2 = seq(from = (lam[2] +1), to = (nrow(data2run)-1), by = 1)/(nrow(data2run)-1)
  # we only consider those that are greater than or equal to 5
  alp1_num = sum(alp_seq1*nrow(data2run) >= 5)
  alp2_num = sum(alp_seq2*nrow(data2run) >= 5)
  expect_equal(nrow(opt_output$all_res),alp1_num+alp2_num)
})

