  drif[0][0] = x2*-3.0+epsilon*(x0*x0)*(1.0/2.0)-epsilon*(x1*x1)*(1.0/2.0)+(epsilon*epsilon)*nu*x2-epsilon*x0*x4-epsilon*x1*x5;
  drif[0][1] = x3*-3.0+(epsilon*epsilon)*nu*x3+epsilon*x0*x1-epsilon*x0*x5+epsilon*x1*x4;
  drif[0][2] = x4*-8.0+(epsilon*epsilon)*nu*x4+epsilon*x0*x2*(3.0/2.0)-epsilon*x1*x3*(3.0/2.0);
  drif[0][3] = x5*-8.0+(epsilon*epsilon)*nu*x5+epsilon*x0*x3*(3.0/2.0)+epsilon*x1*x2*(3.0/2.0);
