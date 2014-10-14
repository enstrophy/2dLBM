#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#ifndef nx
#define nx 256 
#endif
const int ny = nx;

double f[nx][ny][9], f_eq[nx][ny][9];
double u[nx][ny], v[nx][ny], rho[nx][ny];
double uOld[nx][ny], vOld[nx][ny], rhoOld[nx][ny];

const double Re = 1000;
const double L = 1.0;

//const double uTop = 0.1;
//const double tau = 3*uTop*(nx-1)/Re+0.5;
const double tau = 2./3.;
const double uTop = Re*(2.*tau-1.)/(6.*(nx-1.));

const double omega = 1.0/tau;
const long int maxTimeStep = 100000;
const double convTol = 1e-6;
const int printFrequency = 1000;

void dumpTecplot(double [][ny], double [][ny]);
void dumpBinary(double [][ny], double [][ny]);
void dumpCenterline(double [][ny], double [][ny]);
double infinityNorm(double [][ny], double [][ny]);
int checkNan(double [][ny]);

int main() {
  //Setting the parameters

  const double w1 = 4.0/9.0, w2 = 1.0/9.0, w3 = 1.0/36.0;

  std::cout << "Re : " << Re << "\n";
  std::cout << "N : " << nx << "\n";
  std::cout << "uTop : " << uTop << "\n";
  std::cout << "Omega : " << omega << "\n";
  const double Ma = uTop*std::sqrt(3.0);
  std::cout << "Mach number : " <<  Ma << "\n"; //Asisuming dx=dt or c=1 
  if (Ma > 0.1) {
    std::cout << "WARNING : Ma large, Consider using more grid points \n";
  }
  std::cout << std::setprecision(6);
  std::cout << std::scientific;
  
  // Setting initial conditions u=0, v=0, f_i=w_i in the interior
  // and u = uTop on the top wall
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny-1; ++j) {
      f[i][j][0] = w1;
      f[i][j][1] = f[i][j][2] = f[i][j][3] = f[i][j][4] = w2;
      f[i][j][5] = f[i][j][6] = f[i][j][7] = f[i][j][8] = w3;
    }
    // Setting for top wall
    f[i][ny-1][0] = w1;
    f[i][ny-1][1] = f[i][ny-1][2] = f[i][ny-1][3] = f[i][ny-1][4] = w2;
    f[i][ny-1][5] = f[i][ny-1][6] = w3;
    f[i][ny-1][7] = w3-uTop/2.0;
    f[i][ny-1][8] = w3+uTop/2.0;
  }

  // Initializing old varbles to check convergence
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny-1; ++j) {
      rhoOld[i][j] = 1.0;
      uOld[i][j] = 0.0;
      vOld[i][j] = 0.0;
    }
    rhoOld[i][ny-1] = 1.0;
    uOld[i][ny-1] = uTop;
    vOld[i][ny-1] = 0.0;
  }

  for (long int t = 0; t < maxTimeStep; t++) {
    // Propagation Step
    // Updating f(1)
    for (int i = nx-1; i > 0; i--) {
      for (int j = 0; j < ny; j++) {
        f[i][j][1] = f[i-1][j][1]; 
      }
    } 
    // Updating f(2)
    for (int i = 0; i < nx; i++) {
      for (int j = ny-1; j > 0; j--) {
        f[i][j][2] = f[i][j-1][2]; 
      }
    } 
    // Updating f(3)
    for (int i = 0; i < nx-1; i++) {
      for (int j = 0; j < ny; j++) {
        f[i][j][3] = f[i+1][j][3]; 
      }
    } 
    // Updating f(4)
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny-1; j++) {
        f[i][j][4] = f[i][j+1][4]; 
      }
    } 
    // Updating f(5)
    for (int i = nx-1; i > 0; i--) {
      for (int j = ny-1; j > 0; j--) {
        f[i][j][5] = f[i-1][j-1][5]; 
      }
    } 
    // Updating f(6)
    for (int i = 0; i < nx-1; i++) {
      for (int j = ny-1; j > 0; j--) {
        f[i][j][6] = f[i+1][j-1][6]; 
      }
    } 
    // Updating f(7)
    for (int i = 0; i < nx-1; i++) {
      for (int j = 0; j < ny-1; j++) {
        f[i][j][7] = f[i+1][j+1][7]; 
      }
    } 
    // Updating f(8)
    for (int i = nx-1; i > 0; i--) {
      for (int j = 0; j < ny-1; j++) {
        f[i][j][8] = f[i-1][j+1][8]; 
      }
    } 
    // Boundary Conditions
    for (int i = 0; i < nx; i++) {
      // Bottom wall :  Bounce Back
      f[i][0][5] = f[i][0][7];    // f(5) = f(7)
      f[i][0][2] = f[i][0][4];    // f(2) = f(4)
      f[i][0][6] = f[i][0][8];    // f(6) = f(8)
      // Top wall :  Zho-Hu velocity
      double rhoN = f[i][ny-1][0]+f[i][ny-1][1]+f[i][ny-1][3]+2.0*(\
          f[i][ny-1][5]+f[i][ny-1][2]+f[i][ny-1][6]);
      f[i][ny-1][4] = f[i][ny-1][2];
      f[i][ny-1][7] = f[i][ny-1][5]+0.5*(f[i][ny-1][1]-f[i][ny-1][3])-uTop*rhoN/2.0;
      f[i][ny-1][8] = f[i][ny-1][6]-0.5*(f[i][ny-1][1]-f[i][ny-1][3])+uTop*rhoN/2.0;
    }
    for (int j = 1; j < ny-1; j++) {
      // Left Wall : Bounce Back
      f[0][j][5] = f[0][j][7];
      f[0][j][1] = f[0][j][3];
      f[0][j][8] = f[0][j][6];
      // Right Wall : Bounce Back
      f[nx-1][j][7] = f[nx-1][j][5];
      f[nx-1][j][3] = f[nx-1][j][1];
      f[nx-1][j][6] = f[nx-1][j][8];
    }
    // Collision
    // Calculating the physical variables, f_eq, f(t+dt)
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        // Density
        rho[i][j] = 0.0;
        for (int m = 0; m < 9; m++) {
          rho[i][j] += f[i][j][m];
        }
        // Velocity
        u[i][j] = ((f[i][j][1]+f[i][j][5]+f[i][j][8])-\
            (f[i][j][3]+f[i][j][6]+f[i][j][7]))/rho[i][j];
        v[i][j] = ((f[i][j][2]+f[i][j][5]+f[i][j][6])-\
            (f[i][j][4]+f[i][j][7]+f[i][j][8]))/rho[i][j];
        //Equilibrium variables
        double factorVel = 1.5*(u[i][j]*u[i][j]+v[i][j]*v[i][j]);

        f_eq[i][j][1] = w2*rho[i][j]*(1.0+3.0*u[i][j]+4.5*u[i][j]*u[i][j]-\
            factorVel);
        f_eq[i][j][2] = w2*rho[i][j]*(1.0+3.0*v[i][j]+4.5*v[i][j]*v[i][j]-\
            factorVel);
        f_eq[i][j][3] = w2*rho[i][j]*(1.0-3.0*u[i][j]+4.5*u[i][j]*u[i][j]-\
            factorVel);
        f_eq[i][j][4] = w2*rho[i][j]*(1.0-3.0*v[i][j]+4.5*v[i][j]*v[i][j]-\
            factorVel);

        f_eq[i][j][5] = w3*rho[i][j]*(1.0+3.0*(u[i][j]+v[i][j])+4.5*(\
              u[i][j]+v[i][j])*(u[i][j]+v[i][j])-factorVel);
        f_eq[i][j][6] = w3*rho[i][j]*(1.0+3.0*(-u[i][j]+v[i][j])+4.5*(\
              -u[i][j]+v[i][j])*(-u[i][j]+v[i][j])-factorVel);
        f_eq[i][j][7] = w3*rho[i][j]*(1.0-3.0*(u[i][j]+v[i][j])+4.5*(\
              u[i][j]+v[i][j])*(u[i][j]+v[i][j])-factorVel);
        f_eq[i][j][8] = w3*rho[i][j]*(1.0+3.0*(u[i][j]-v[i][j])+4.5*(\
              u[i][j]-v[i][j])*(u[i][j]-v[i][j])-factorVel);

        f_eq[i][j][0] = w1*rho[i][j]*(1.0-factorVel);
        
        // Updating f
        for (int m = 0; m < 9; m++) {
          f[i][j][m] = omega*f_eq[i][j][m]+(1.0-omega)*f[i][j][m];
        }
      }
    }
    //Check for Nan
    if (checkNan(u) || checkNan(v) || checkNan(rho)) {
      std::cout << "Iteration : " << t+1 << " ...exiting\n";
      return 1;
    }
    //Check for convergence
    double rhoConv = infinityNorm(rho, rhoOld);
    double uConv = infinityNorm(u, uOld);
    double vConv = infinityNorm(v, vOld);
    // std::cout << "\n";
    if ( (rhoConv < convTol) && (uConv < convTol) && (vConv < convTol) ) {
      std::cout << "\nConverged in " << t+1 << " iterations \n";
      break;
    }
    // Update old variables
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        rhoOld[i][j] = rho[i][j];
        uOld[i][j] = u[i][j];
        vOld[i][j] = v[i][j];
      }
    }
    if (!((t+1)%printFrequency)) {
      std::cout << "Iteration : " << t+1 << "\t";
      std::cout << rhoConv << "\t" << uConv << "\t" << vConv << "\n";
      // int index = nx/2-1;
      // std::cout << rho[index][index]<< "\t" << u[index][index] << "\t"
      //          << v[index][index] << "\n";
    }
  }
  dumpTecplot(u, v);
  dumpBinary(u, v);
  dumpCenterline(u, v);
  return 0;
}

double infinityNorm(double array1[][ny], double array2[][ny]) {
  double norm = 0.0;
  int ii = 0, jj = 0;
  norm = std::fabs(array1[0][0]-array2[0][0]);
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      if (norm < std::fabs(array1[i][j]-array2[i][j])) {
        norm = std::fabs(array1[i][j]-array2[i][j]);
        ii = i; jj = j;
      }
    }
  }
  // std::cout << ii << "\t" << jj << "\t" << array1[ii][jj] << "\t" << array2[ii][jj] << "\n";
  return norm;
}
void dumpTecplot(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1.open("lidCavityASCII.dat");
  f1 << std::scientific;
  f1.precision(4);
  f1 << "TITLE = Re=" << Re << "1000 Lid Driven Cavity" << std::endl \
     << "VARIABLES = X, Y, U, V" << std::endl;
  f1 << "Zone I =" << nx << " J=" << ny << "F=POINT" << std::endl;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double x = (i+0.5)/nx;
      double y = (j+0.5)/ny;
      f1 << x << "\t" << y << "\t" << u[i][j]/uTop << "\t" << v[i][j]/uTop << "\n";
    }
  }
  f1.close();
}
void dumpBinary(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1.open("lidCavityResult.dat", std::ios::binary);
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double x = (i+0.5)/nx;
      double y = (j+0.5)/ny;
      double vmag = std::sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
      vmag /= uTop;
      f1.write((char*)&x, sizeof(double));
      f1.write((char*)&y, sizeof(double));
      f1.write((char*)&vmag, sizeof(double));
    }
  }
  f1.close();
}
void dumpCenterline(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1 << std::setprecision(6);
  f1 << std::scientific;
  f1.open("ldcVerticalCenterline.dat");
  int index = int(nx/2);
  if ((nx%2) != 0) {
    for (int i = 0; i < ny; ++i) {
      f1 << i/double(ny-1) << "\t" << u[index][i]/uTop << "\n";
    }
  } else {
    for (int i = 0; i < ny; ++i) {
      f1 << i/double(ny-1) << "\t" << (u[index][i]+u[index-1][i])/(2.0*uTop) << "\n";
    }
  }
  f1.close();
  f1.open("ldcHorizontalCenterline.dat");
  index = int(ny/2);
  if ((ny%2) != 0) {
    for (int i = 0; i < nx; ++i) {
      f1 << i/double(nx-1) << "\t" << v[i][index]/uTop << "\n";
    }
  } else {
    for (int i = 0; i < nx; ++i) {
      f1 << i/double(nx-1) << "\t" << (v[i][index]+v[i][index-1])/(2.0*uTop) << "\n";
    }
  }
  f1.close();
}
int checkNan(double u[][ny]) {
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      if (u[i][j] != u[i][j]) {
        std::cout << "Nan found at " << i << "\t" << j << "\n";
        return 1;
      }
    }
  }
  return 0;
}
