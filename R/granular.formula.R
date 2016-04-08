#' @export
atm = list(
  e = list(m = 0.000548579909),
  p = list(m = 1.00727647),
  n = list(m = 1.008671087),
  
  n14 = list(m = 14.0030740048, a = 0.99632),
  n15 = list(m = 15.0001088982, a = 0.003268),
  
  c12 = list(m = 12.0000, a = 0.9893),
  c13 = list(m = 13.0033548378, a = 0.0107),
  c14 = list(m = 14.003241989, a = 0),
  
  h1 = list(m = 1.00782503207, a = 0.999885),
  h2 = list(m = 2.0141017778, a = 0.000115),
  
  o16 = list(m = 15.994915, a = 0.99757),
  o17 = list(m = 16.999132, a = 0.00038),
  o18 = list(m = 17.999160, a = 0.00205),
  
  li6 = list(m = 6.015122, a = 0.0759),
  li7 = list(m = 7.016004, a = 0.9241),
  
  na23 = list(m = 22.989770, a = 1),
  
  k39 = list(m = 38.963707, a = 0.932581),
  k41 = list(m = 40.961826, a = 0.067302),
  
  cl35 = list(m = 34.968853, a = 0.7578),
  cl37 = list(m = 36.965903, a = 0.2422),
  
  br79 = list(m = 78.918338, a = 0.5069),
  br81 = list(m = 80.916291, a = 0.4931),
  
  mg24 = list(m = 23.985042, a = 0.7899),
  mg25 = list(m = 24.985837, a = 0.1000),
  mg26 = list(m = 25.982593, a = 0.1101),
  
  s32 = list(m = 31.972071, a = 0.9493),
  s33 = list(m = 32.971458, a = 0.0076),
  s34 = list(m = 33.967867, a = 0.0429),
  
  f19 = list(m = 18.998403, a = 1),
  
  p31 = list(m = 30.973762, a = 1),
  
  si28 = list(m = 27.976927, a = 0.922297),
  si29 = list(m = 28.976495, a = 0.046832),
  si30 = list(m = 29.973770, a = 0.030872)
)


#' @export
M.iso = as.data.frame(matrix(c(
  atm$c13$m - atm$c12$m, atm$c13$a, 0.071, 0, 1,
  atm$n15$m - atm$n14$m, atm$n15$a, 0.040, 0, 1,
  atm$o18$m - atm$o16$m, atm$o18$a, 0.040, 0, 1,
  
  atm$s33$m - atm$s32$m, atm$s33$a, 0.020, 0, 1,
  atm$s34$m - atm$s32$m, atm$s34$a, 0.020, 0, 1,
  
  atm$cl37$m - atm$cl35$m, atm$cl37$a, 0.020, 0, 1,
  atm$br81$m - atm$br79$m, atm$br81$a, 0.008, 0, 1,
  
  atm$si29$m - atm$si28$m, atm$si29$a, 0.008, 0, 1,
  atm$si30$m - atm$si28$m, atm$si30$a, 0.008, 0, 1,
  
  atm$k41$m - atm$k39$m, atm$k41$a, 0.008, 0, 1
),
ncol = 5, byrow=T,
dimnames = list(
  c("C12-13", "N14-15", "O16-18", "S32-33", "S32-34", "Cl35-37", "Br79-81", "Si28-29", "Si28-30", "K41-39"), 
  c("m", "ab", "npm", "z", "d")
)
))


# Name, charge added, mass of adduct
#' @export
M.z = as.data.frame(matrix(c(
  +1, atm$h1$m, 0.9, 0,
  +1, atm$na23$m, 0.8, 0,
  +1, atm$k39$m, 0.7, 0,
  -1, atm$cl35$m, 0.5, 0,
  -1, atm$br79$m, 0.2, 0
),
ncol = 4, byrow=T,
dimnames = list(
  c(
    "H+", 
    "Na+", 
    "K+", 
    "Cl-",
    "Br-"
  ), 
  c("z", "m", "p", "d")
)
))

M.z.specs = list(
  c(1,100,2,0.016),
  c(23,100),
  c(39,100,40,0.0129,41,7.221),
  c(35,100,37,32.3977),
  c(79,100,81,97.5114)
)
M.z.specs = lapply(M.z.specs, function(x) { data.frame(matrix(x, ncol=2, byrow=T, dimnames = list(NULL, c("m", "i")))) %>% {.$i = .$i / max(.$i); .} })
names(M.z.specs) = M.z$id


#' @export
M.n = as.data.frame(matrix(c(
  0, -(atm$h1$m*2+atm$o16$m), 0.6, 1,
  0, -(atm$c12$m + atm$o16$m*2), 0.4, 1,
  0, -(atm$n14$m + atm$h1$m*3), 0.4, 1,
  0, +(atm$c12$m + atm$o16$m*2 + atm$h1$m*2), 0.6, 1,
  0, +(atm$c12$m*2 + atm$h1$m*4 + atm$o16$m*2), 0.6, 1,
  0, +(atm$c12$m*2 + atm$h1$m*1 + atm$o16$m*2 + atm$f19$m*3), 0.1, 1,
  0, +(atm$c12$m*2 + atm$h1$m*3 + atm$n14$m), 0.3, 1,
  0, +(atm$c12$m + atm$h1$m*4 + atm$o16$m), 0.3, 1,
  0, -(atm$c12$m + atm$o16$m), 0.3, 1,
  0, +(atm$h1$m*3 + atm$p31$m + atm$o16$m*4), 0.1, 1,
  0, +(atm$si28$m + atm$h1$m*2 + atm$o16$m*3), 0.1, 1,
  0, +(atm$si28$m + atm$h1$m*4 + atm$o16$m*4), 0.1, 1,
  0, +(atm$si28$m + atm$h1$m*6 + atm$o16$m + atm$c12$m*2), 0.1, 1
),
ncol=4, byrow=T,
dimnames= list(
  c("-H2O", 
    "-CO2", 
    "-NH3",
    "+HCOOH",
    "+CH3COOH",
    "+CF3COOH",
    "+CH3CN", 
    "+CH3OH",
    "-CO",
    "+H3PO4",
    "+SiO3H2",
    "+SiO4H4",
    "+SiC2H6O"
  ),
  c("z", "m", "p", "d")
)
))


M.iso$id = 1:(length(M.iso$m))
M.z$id = 1:(length(M.z$m)) + max(M.iso$id)
M.n$id = 1:(length(M.n$m)) + max(M.z$id)



M.n.specs = list(
  c(18,100,19,0.0721,20,0.2005,21,0.0001),
  c(44,100,45,1.1618,46,0.4018,47,0.0045,48,0.0004),
  c(17,100,18,0.4093,19,0.0002),
  c(46,100,47,1.1938,48,0.4022,49,0.0046,50,0.0004),
  c(60,100,61,2.3073,62,0.4158,63,0.0091,64,0.0004),
  c(114,100,115,2.2593,116,0.4147,117,0.0089,118,0.0004),
  c(41,100,42,2.5725,43,0.0207),
  c(32,100,33,1.1857,34,0.2016,35,0.0023),
  c(28,100,29,1.1217,30,0.2009,31,0.0022),
  c(98,100,99,0.2084,100,0.8021,101,0.0013,102,0.0024),
  c(78,100,79,5.2157,80,3.9703,81,0.0362,82,0.0214,83,0.0001),
  c(96,100,97,5.2878,98,4.1746,99,0.0496,100,0.0294,101,0.0001,102,0.0001),
  c(74,100,75,7.3627,76,3.6927,77,0.0927,78,0.0074,79,0.0001)
)
M.n.specs = lapply(M.n.specs, function(x) { data.frame(matrix(x, ncol=2, byrow=T, dimnames = list(NULL, c("m", "i")))) %>% {.$i = .$i / max(.$i); .} })
names(M.n.specs) = M.n$id