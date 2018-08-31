DAISIE_loglik_rhs = function(t,x,pars)

{

  lx = (length(x) - 1)/2
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  kk = pars[6]
  ddep = pars[7]

  nn = -2:(lx+2*kk+1)
  lnn = length(nn)
  nn = pmax(rep(0,lnn),nn)

  xx1 = c(0,0,x[1:lx],0)
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 = x[2 * lx + 1]

  nil2lx = 3:(lx + 2)

  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2

  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk

  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

  dx1 = laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
        lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
        -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
        -gamvec[il3] * xx1[ix3]

  dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
  dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)

  dx2 = gamvec[il3] * xx1[ix3] +
      lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
      -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
      -laavec[il3 + 1] * xx2[ix3]

  dx3 = -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] + muvec[il3[1]]) * xx3

  return(list(c(dx1,dx2,dx3)))
}
