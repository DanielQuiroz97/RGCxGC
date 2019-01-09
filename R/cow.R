setGeneric('histc',
           def = function(values, edges){
             standardGeneric('histc')
           })
setMethod('histc',
          signature = 'numeric',
          definition =  function(values, edges) {
            ledges <- length(edges)
            nlow <- length(values[values < edges[1]])
            nhigh <- length(values[values > edges[ledges]])
            bin <- NULL
            for (i in seq_along(edges)[-ledges]) {
              bin_values <- values[values >= edges[i] & values < edges[i + 1]]
              position <- match(bin_values, values)
              bin[position] <- i
            }
            bin <- bin[!is.na(bin)]
            upper <- match(edges[ledges], values)
            if (!is.na(upper)) {
              bin[length(bin) + 1] <- i + 1
            }
            if (nlow > 0) {
              bin <- c(rep(0, nlow), bin)
            }
            if (nhigh > 0) {
              bin <- c(bin, rep(0, nlow))
            }
            return(bin)
          })

setGeneric(name = "InterpCoeff",
           def = function(n, nprime, offs, rtn){
             standardGeneric("InterpCoeff")
           })
setMethod(f = "InterpCoeff", signature = "numeric",
          definition = function(n, nprime, offs, rtn) {
            p <- length(nprime)
            q <- n - 1
            coeff <- matrix(nrow = p, ncol = n)
            index <- matrix(nrow = p, ncol = n)
            for (i in seq(1, p)) {
              pp <- seq(1, nprime[i])
              p <- seq(0, q) * (nprime[i] - 1) / q + 1
              k <- histc(p, pp)
              k[p >= nprime[i]] <- nprime[i] - 1
              coeff[i, ] <- (p - pp[k])
              index[i, ] <- k - offs[i]
            }
            switch(rtn, coeff = return(coeff), index = return(index))
          })

setGeneric(name = "cow",
           def = function(Ta, X, seg, Slack, Options = c(0, 1, 0, 0, 0)){
             standardGeneric('cow')
           })
setMethod(f = "cow", signature = "numeric",
          definition = function(Ta, X, seg, Slack, Options = c(0, 1, 0, 0, 0)) {
            # Initialise
            if (is.matrix(X)) {
              dim_x <- length(X)
            } else {
              dim_x <- c(1, length(X))
            }
            np_t <- length(Ta)
            #### Initialise segments ####
            seg <- floor(seg)
            pred_bound <- length(seg) > 1
            if (pred_bound) {
              if (!all(seg[, 1] == 1 & seg[length(seg)] == length(Ta))) {
                stop("End points must be equal to 1 and to
                     the length of the target")
              }
              len_segs <- apply(seg, 1, diff)
              if (!all(len_segs > 2)) {
                stop("segments must contain at least two points")
              }
              len_segs <- t(len_segs)
              n_seg <- ncol(len_segs)
            } else {
              if (seg > min(c(dim_x[2], np_t))) {
                stop("segment length is larger than length of the signal")
              }
              if (Options[[3]]) {
                n_seg <- floor( (np_t - 1) / seg)
                len_segs <- matrix( floor( (np_t - 1 ) / n_seg), nrow = 1)
                len_segs[2, ] <- floor( (dim_x[2] - 1) / n_seg)
                print("segment lengh adjusted to the best cover the remainders")
              } else {
                n_seg <- floor( (np_t - 1 ) / (seg - 1))
                tmp_segs <- rep(seg - 1, n_seg)
                len_segs <- matrix(c(tmp_segs, tmp_segs), nrow = 2)
                if (floor( ( dim_x[2] - 1) / (seg - 1)) != n_seg) {
                  stop("For non-fixed segment lengths the target and teh signal
             do not have the same number o fsegments. Try option 3 set to T")
                }
              }
              tmp <- (np_t - 1) %% len_segs[1, 1]
              if (tmp > 0) {
                len_segs[1, n_seg] <- len_segs[1, n_seg] + tmp
                if (Options[[1]]) {
                  print(paste0("segments: ", len_segs[1, 1] + 1,
                               " Points: ", n_seg - 1))
                }
              } else {
                if (Options[[1]]) {
                  print(paste0("segments: ", len_segs[2, 1] + 1,
                               " Points x: ", n_seg))
                }
              }
              tmp <- (dim_x[2] - 1) %% len_segs[2, 1]
              if (tmp > 0) {
                len_segs[2, n_seg] <- len_segs[2, 1] + tmp
              }
            }
            bt <- cumsum(c(1, len_segs[1, ]))
            bp <- cumsum(c(1, len_segs[2, ]))
            Warping <- matrix(nrow = dim_x[1], ncol = (n_seg + 1))
            #### Chech Slack ####
            if (length(Slack) > 1) {
              if (length(Slack) <= n_seg) {
                stop("The number of slack parameters is not equal to
           the number of optimised segments")
              }
              stop("Multiple slacks have not been implemented yet")
            }
            Slacks_vec <- seq(-Slack, Slack)
            #### Set feasible points for boundaies ####
            Bounds <- matrix(rep(1, 2 * (n_seg + 1)), nrow = 2)
            # Slope constrints
            offs_tmp <- Slack * seq(0, n_seg)
            offs <- matrix(c(-offs_tmp, offs_tmp), nrow = 2, byrow = T)
            Bounds_ind <- seq(1, n_seg + 1)
            Bounds_a2 <- matrix( rep( bp[Bounds_ind], 2), nrow = 2, byrow = T)
            Bounds_a <- Bounds_a2 + offs
            offs_tmpb <- matrix(c(-rev(offs_tmp), rev(offs_tmp)),
                                nrow = 2, byrow = T)
            Bounds_b <- Bounds_a2 + offs_tmpb
            Bounds[1, ] <- apply(matrix(c(Bounds_a[1, ], Bounds_b[1, ]),
                                        nrow = 2, byrow = T), 2, max)
            Bounds[2, ] <- apply(matrix(c(Bounds_a[2, ], Bounds_b[2, ]),
                                        nrow = 2, byrow = T), 2, min)
            #### Calculate de first derivate for interpolation ####
            Xdiff <- diff(X)
            #### Calculate coefficients and indexes for interpolation ####
            Int_Coeff <- vector("list", n_seg)
            Int_Index <- Int_Coeff
            if (!pred_bound) {
              for (i in seq(1, n_seg - 1)) {
                nprm <- len_segs[2, 1] + Slacks_vec + 1
                Int_Coeff[[i]] <- InterpCoeff(n = len_segs[1, 1] + 1,
                                              nprime = nprm,
                                              offs = Slacks_vec,
                                              rtn = "coeff")
                Int_Index[[i]] <- InterpCoeff(n = len_segs[1, 1] + 1,
                                              nprime = nprm,
                                              offs = Slacks_vec, rtn = "index")
              }
              nprm2 <-  len_segs[2, n_seg] + Slacks_vec + 1
              Int_Coeff[[n_seg]] <- InterpCoeff(n = len_segs[1, n_seg] + 1,
                                                nprime = nprm2,
                                                offs = Slacks_vec,
                                                rtn = "coeff")
              Int_Index[[n_seg]] <- InterpCoeff(n = len_segs[1, n_seg] + 1,
                                                nprime = nprm2,
                                                offs = Slacks_vec,
                                                rtn = "index")
            } else {
              for (i in seq(1, n_seg)) {
                nprm <- len_segs[2, i] + Slacks_vec + 1
                Int_Coeff[[i]] <- InterpCoeff(n = len_segs[1, i] + 1,
                                              nprime = nprm,
                                              offs = Slacks_vec,
                                              rtn = "coeff")
                Int_Index[[i]] <- InterpCoeff(n = len_segs[1, i] + 1,
                                              nprime = nprm,
                                              offs = Slacks_vec, rtn = "index")
              }
            }
            #### Dynamic Programming Section ####
            table_index <- cumsum(c(0, Bounds[2, ] - Bounds[1, ] + 1))
            Table <- matrix(0, nrow = 3, ncol = table_index[n_seg + 2])
            Table[2, 2:ncol(Table)] <- -Inf
            for (i in seq(1, n_seg + 1)) {
              v <- seq(Bounds[1, i], Bounds[2, i])
              Table[1, seq(table_index[i] + 1, table_index[i + 1])] <- v
            }
            # Forward phase
            for (i in seq(1, n_seg)) {
              a <- Slacks_vec + len_segs[2, i]
              b <- table_index[i] + 1 - Bounds[1, i]
              c <- len_segs[1, i] + 1
              counting <- 1
              node_z <- table_index[i + 2]
              node_a <- table_index[i + 1] + 1
              bound_k_table <- matrix(nrow = 2, ncol = node_z - node_a + 1)
              int_index_seg <- t(Int_Index[[i]]) - (len_segs[2, i] + 1)
              int_coeff_seg <- t(Int_Coeff[[i]])
              Tseg <- Ta[seq(bt[i], bt[i + 1])]
              Tseg_centered <- Tseg - sum(Tseg) / length(Tseg)
              Norm_Tseg_cen <- norm(Tseg_centered, type = "2")
              for (j in seq(node_a, node_z)) {
                prec_nodes <- Table[1, j] - a
                allowed_arcs <- prec_nodes >= Bounds[1, i] & 
                  prec_nodes <= Bounds[2, i]
                nodes_tablepointer <- b + prec_nodes[allowed_arcs]
                n_aa <- sum(allowed_arcs)
                if (n_aa != 0) {
                  Index_Node <- Table[1, j] + int_index_seg[, allowed_arcs]
                  coeff_b <- sapply(int_coeff_seg[, allowed_arcs],
                                    function(x) x)
                  Xi_seg <- X[Index_Node]
                  Xi_diff <- Xdiff[Index_Node]
                  toreshape <- matrix(c(coeff_b, Xi_diff), nrow = 2, byrow = T)
                  toreshape <- apply(toreshape, 2, prod) + Xi_seg
                  Xi_seg <- matrix(toreshape, nrow = c, ncol = n_aa * dim_x[1])
                  Xi_seg_mean <- colSums(Xi_seg) / nrow(Xi_seg)
                  Norm_Xi_seg_cen <- sqrt(colSums(Xi_seg ^ 2) - nrow(Xi_seg) *
                                            Xi_seg_mean ^ 2)
                  CCs_Node <- (as.numeric(Tseg_centered) %*% Xi_seg) /
                    (Norm_Tseg_cen %*% Norm_Xi_seg_cen)
                  CCs_Node <- ifelse(is.finite(CCs_Node), CCs_Node, 0)
                  CCs_Node <- matrix(CCs_Node, nrow = n_aa, ncol = dim_x[1])
                  if (Options[[2]]) {
                    Cost_Fun <- matrix(Table[2, nodes_tablepointer],
                                       nrow = n_aa) + CCs_Node
                  } else {
                    Cost_Fun <- matrix(Table[2, nodes_tablepointer],
                                       nrow = n_aa) + CCs_Node ^ Options[2]
                  }
                  ind <- max(Cost_Fun)
                  pos <- match(ind, Cost_Fun)
                  bound_k_table[1, counting] <- ind
                  bound_k_table[2, counting] <- nodes_tablepointer[pos]
                  counting <- counting + 1
                }
              }
              Table[2:3, seq(node_a, node_z)] <- bound_k_table
            }
            for (i in seq(1, dim_x[1])) {
              Pointer <- ncol(Table)
              Warping[i, n_seg + 1] <- dim_x[2]
              for (j in seq(n_seg, 1, -1)) {
                Pointer <- Table[3, Pointer]
                Warping[i, j] <- Table[1, Pointer]
              }
            }
            Xwarped <- NULL
            for (i in seq(1, n_seg)) {
              indt <- seq(bt[i], bt[i + 1])
              lent <- bt[i + 1] - bt[i]
              for (j in seq(1, dim_x[1])) {
                ind_x <- seq(Warping[j, i], Warping[j, i + 1])
                len_x <- Warping[j, i + 1] - Warping[j, i]
                Xwarped[indt] <- approx(x = ind_x - Warping[j, i] + 1,
                                        y = X[ind_x],
                                        xout = seq(0,
                                                   lent) / lent * len_x + 1)$y
              }
            }
            return(list(Warping = Warping, XWarped = Xwarped))
          })