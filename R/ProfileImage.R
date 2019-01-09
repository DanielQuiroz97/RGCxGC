setGeneric(name = "ProfileImage", 
           def = function(in_image, in_pdim, in_t, in_width){
             standardGeneric("ProfileImage")
           })
setMethod(f = "ProfileImage",
          signature = "matrix",
          definition = function(in_image, in_pdim, in_t, in_width) {
            row <- nrow(in_image)
            col <- ncol(in_image)
            if (in_pdim == 1) {
              fromseq <- in_t - in_width
              fromseq <- ifelse(fromseq < 1, 1, fromseq)
              toseq <- in_t + in_width
              toseq <- ifelse(toseq > row, row, toseq)
              tmp_idx <- seq(fromseq, toseq)
              tmp_x <- in_image[tmp_idx, ]
              tmp_w <- 0.75 * (1 - ((tmp_idx - in_t) / in_width) ^ 2)
              tensor <- kronecker(matrix(1, ncol = col), tmp_w)
              tmp_tensor <- tmp_x * tensor
              ret_profile <- colSums(tmp_tensor)
            } else {
              fromseq <- in_t - in_width
              fromseq <- ifelse(fromseq < 1, 1, fromseq)
              toseq <- in_t + in_width
              toseq <- ifelse(toseq > col, col, toseq)
              tmp_idx <- seq(fromseq, toseq)
              tmp_x <- in_image[, tmp_idx]
              tmp_w <- 0.75 * (1 - ((tmp_idx - in_t) / in_width) ^ 2)
              tensor <- kronecker(matrix(1, nrow = row), tmp_w)
              tensor <- matrix(tensor, nrow = row, byrow = T)
              tmp_tensor <- tmp_x * tensor
              ret_profile <- rowSums(tmp_tensor)
            }
            return(ret_profile)
          })