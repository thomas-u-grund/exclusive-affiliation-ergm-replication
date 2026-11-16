# Term sets of the football models (same as nwergm: stage 1, 3m, 3ms2). Evaluate with `fb` loaded (ergm_excl_data.R).
excl_formula <- function(model, fb) {
  e <- list2env(fb, parent = asNamespace("ergm.exclusive")); e$ergm <- NULL
  .sz <- 10 * exp(fb$wsize); ag2 <- fb$isag == 1 & .sz >= 2 - 1e-6; .mu <- mean(log(.sz[ag2])); .sg <- sd(log(.sz[ag2]))
  e$zsize <- ifelse(fb$isag == 1, (log(.sz) - .mu) / .sg, 0); e$hinge <- ifelse(fb$isag == 1, pmax(0, log(.sz) - log(10)) / .sg, 0)
  e$sml <- fb$small * (.sz < 10 - 1e-6); e$mid <- fb$small * (.sz >= 10 - 1e-6)
  .c1 <- as.numeric(Sys.getenv("CUT1", "10")); .c2 <- as.numeric(Sys.getenv("CUT2", "30"))   # size-class cut-offs of model 3s3c (sensitivity)
  e$smlc <- fb$isag * (.sz < .c1 - 1e-6); e$midc <- fb$isag * (.sz >= .c1 - 1e-6 & .sz < .c2 - 1e-6); e$lrgc <- fb$isag * (.sz >= .c2 - 1e-6)
  rhs <- switch(model,
    "1"   = "stack() + homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country')",
    "3m"  = paste("stack() + homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                  "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                  "edgecov(pc_home) + edgecov(pc_val)"),
    "3ms2" = paste("homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                   "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                   "edgecov(pc_home) + edgecov(pc_val) + coloc(small, club_id, 'small') + coloc(large, club_id, 'large')"),
    "3sz" = paste("stack() + coloc(wsize, club_id, 'size') + homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                  "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                  "edgecov(pc_home) + edgecov(pc_val)"),
    "3s3" = paste("homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                  "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                  "edgecov(pc_home) + edgecov(pc_val) + coloc(sml, club_id, 'under10') + coloc(mid, club_id, '10to29') + coloc(large, club_id, 'large')"),
    "3s3c" = paste("homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                  "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                  "edgecov(pc_home) + edgecov(pc_val) + coloc(smlc, club_id, 'under10') + coloc(midc, club_id, '10to29') + coloc(lrgc, club_id, 'large')"),
    "3d" = paste("homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                 "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                 "edgecov(pc_home) + edgecov(pc_val) + coloc(small, club_id, 'small') + coloc(large, club_id, 'large') + dstack(small, 'small')"),
    "3szs" = paste("stack() + coloc(zsize, club_id, 'zsize') + homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                   "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                   "edgecov(pc_home) + edgecov(pc_val)"),
    "3spl" = paste("stack() + coloc(zsize, club_id, 'zsize') + coloc(hinge, club_id, 'hinge') + homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                   "acov(ac_home, 'home') + acov(ac_val, 'val') + acov(ac_mkt, 'mkt') + edgecov(ap_corr) + edgecov(ap_val) +",
                   "edgecov(pc_home) + edgecov(pc_val)"),
    "3s3nm" = paste("homcentre(isag, nat, 'agency') + homcentre(iscl, nat, 'club') + coloc(isag, key_cty, 'country') +",
                  "acov(ac_home, 'home') + acov(ac_val, 'val') + edgecov(ap_corr) + edgecov(ap_val) +",
                  "edgecov(pc_home) + edgecov(pc_val) + coloc(sml, club_id, 'under10') + coloc(mid, club_id, '10to29') + coloc(large, club_id, 'large')"),
    stop("unknown model"))
  f <- as.formula(paste("fb_net ~", rhs)); e$fb_net <- fb$net; environment(f) <- e; f
}
