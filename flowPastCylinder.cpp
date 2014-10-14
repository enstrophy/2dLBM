#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifndef nx
#define nx 400
#endif

#ifndef ny
#define ny 100
#endif

#ifndef maxTimeStep 
#define maxTimeStep 400000
#endif

const int circleX = nx/5;
const int circleY = ny/2+2;  // Mounted asymmetrically to start shedding
const int circleR = ny/10+1;
// const int circleR = 0;
const int circleR2 = circleR*circleR;

const double uMax = 0.1;
const double Re = 100;
const double nu = uMax*2.*circleR/Re;
const double tau = (3.*nu+0.5);
const double omega = 1./tau;

const int printFrequency = 1000;

// D3Q9 Constants
const double weights[] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., \
  1./36., 1./36., 1./36.};
const double cx[] = {0., 1., 0., -1., 0., 1., -1., -1., 1.};
const double cy[] = {0., 0., 1., 0., -1., 1., 1., -1., -1.};
const int opp[] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

double f_in[nx][ny][9], f_eq[nx][ny][9], f_out[nx][ny][9];
double u[nx][ny], v[nx][ny], rho[nx][ny];

void dumpTecplot(double [][ny], double [][ny]);
void dumpTecplotMag(double [][ny], double [][ny]);
void dumpBinary(double [][ny], double [][ny]);
int checkNan(double [][ny]);

int main() {
  std::cout << "Re : " << Re << "\n";
  std::cout << "Nx : " << nx << "\n";
  std::cout << "Ny : " << ny << "\n";
  std::cout << "Omega : " << omega << "\n";
  const double Ma = uMax*std::sqrt(3.0);
  std::cout << "Mach number : " <<  Ma << "\n"; //Asisuming dx=dt or c=1 
  if (Ma > 0.1) {
    std::cout << "WARNING : Ma large, Consider using more grid points \n";
  }
  std::cout << std::setprecision(6);
  std::cout << std::scientific;

  // Initial conditions
  const double width = ny-2.;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      double vl = 0.;
      double y = j-0.5;
      double ul = uMax*4.*(y*width-y*y)/(width*width);
      double rhol = 1.;
      double uMag = 3.*(ul*ul+vl*vl)/2.;
      for (int m =0; m < 9; m++) {
        double t2 = 3.*(cx[m]*ul+cy[m]*vl);
        f_in[i][j][m] = rhol*weights[m]*(1.+t2+t2*t2/2.-uMag);
      }
    }
  }
  // Main loop, time iteration
  for (int t = 0; t < maxTimeStep; ++t) {
    // Calculating macroscopic variables
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        // Density
        rho[i][j] = 0.0;
        for (int m = 0; m < 9; m++) {
          rho[i][j] += f_in[i][j][m];
        }
        // Velocity
        u[i][j] = ((f_in[i][j][1]+f_in[i][j][5]+f_in[i][j][8])-\
            (f_in[i][j][3]+f_in[i][j][6]+f_in[i][j][7]))/rho[i][j];
        v[i][j] = ((f_in[i][j][2]+f_in[i][j][5]+f_in[i][j][6])-\
            (f_in[i][j][4]+f_in[i][j][7]+f_in[i][j][8]))/rho[i][j];
      }
    }
    // Dirichlet boundary conditions
    // Inlet : Poiseuille profile; Zou-He Velocity Inlet
    for (int j = 1; j < ny-1; ++j) {
      double y = j-0.5;
      double vl = 0.;
      double ul = 4.*uMax*(y*width-y*y)/(width*width);
      double rhol = (f_in[0][j][0]+f_in[0][j][2]+f_in[0][j][4]+\
          2.0*(f_in[0][j][3]+f_in[0][j][6]+f_in[0][j][7]))/(1.-ul);
      rho[0][j] = rhol;
      u[0][j] = ul;
      v[0][j] = vl;
      f_in[0][j][1] = f_in[0][j][3]+2.*rhol*ul/3.;
      f_in[0][j][5] = f_in[0][j][7]+0.5*(f_in[0][j][4]-f_in[0][j][2])+\
                      0.5*rhol*vl+rhol*ul/6.;
                      // rhol*ul/6.;
      f_in[0][j][8] = f_in[0][j][6]-0.5*(f_in[0][j][4]-f_in[0][j][2])-\
                      0.5*rhol*vl+rhol*ul/6.;
                      // rhol*ul/6.;
    }
    // Outlet : Zou-He Pressure Outlet
    for (int j = 1; j < ny-1; ++j) {
      double vl = 0.;
      double rhol = 1.0;
      int out = nx-1;
      double ul = -1.0+(f_in[out][j][0]+f_in[out][j][2]+f_in[out][j][4]+\
          2.0*(f_in[out][j][1]+f_in[out][j][5]+f_in[out][j][8]))/rhol;
      rho[out][j] = rhol;
      u[out][j] = ul;
      v[out][j] = vl;
      f_in[out][j][3] = f_in[out][j][1]-2.*rhol*ul/3.;
      f_in[out][j][7] = f_in[out][j][5]-0.5*(f_in[out][j][4]-f_in[out][j][2])-\
                      0.5*rhol*vl-rhol*ul/6.;
                      // rhol*ul/6.;
      f_in[out][j][6] = f_in[out][j][8]+0.5*(f_in[out][j][4]-f_in[out][j][2])+\
                      0.5*rhol*vl-rhol*ul/6.;
                      // -rhol*ul/6.;
    }
    // Collision
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double uMag = 3.*(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.;
        for (int m =0; m < 9; m++) {
          double t2 = 3.*(cx[m]*u[i][j]+cy[m]*v[i][j]);
          f_eq[i][j][m] = rho[i][j]*weights[m]*(1.+t2+t2*t2/2.-uMag);
          f_out[i][j][m] = f_in[i][j][m]-omega*(f_in[i][j][m]-f_eq[i][j][m]);
        }
      }
    }
    // Bounce back at circle
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        int position = std::pow(i-circleX, 2)+std::pow(j-circleY, 2);
        if ((j==0) || (j==ny-1) || (position <= circleR2)) {
          for (int m =0; m < 9; m++) {
            f_out[i][j][m] = f_in[i][j][opp[m]];
          }
        }
      }
    }
    // Streaming
    // Updating f(0)
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        f_in[i][j][0] = f_out[i][j][0]; 
      }
    } 
    // Updating f(1)
    for (int i = nx-1; i > 0; i--) {
      for (int j = 0; j < ny; j++) {
        f_in[i][j][1] = f_out[i-1][j][1]; 
      }
    } 
    // Updating f(2)
    for (int i = 0; i < nx; i++) {
      for (int j = ny-1; j > 0; j--) {
        f_in[i][j][2] = f_out[i][j-1][2]; 
      }
    } 
    // Updating f(3)
    for (int i = 0; i < nx-1; i++) {
      for (int j = 0; j < ny; j++) {
        f_in[i][j][3] = f_out[i+1][j][3]; 
      }
    } 
    // Updating f(4)
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny-1; j++) {
        f_in[i][j][4] = f_out[i][j+1][4]; 
      }
    } 
    // Updating f(5)
    for (int i = nx-1; i > 0; i--) {
      for (int j = ny-1; j > 0; j--) {
        f_in[i][j][5] = f_out[i-1][j-1][5]; 
      }
    } 
    // Updating f(6)
    for (int i = 0; i < nx-1; i++) {
      for (int j = ny-1; j > 0; j--) {
        f_in[i][j][6] = f_out[i+1][j-1][6]; 
      }
    } 
    // Updating f(7)
    for (int i = 0; i < nx-1; i++) {
      for (int j = 0; j < ny-1; j++) {
        f_in[i][j][7] = f_out[i+1][j+1][7]; 
      }
    } 
    // Updating f(8)
    for (int i = nx-1; i > 0; i--) {
      for (int j = 0; j < ny-1; j++) {
        f_in[i][j][8] = f_out[i-1][j+1][8]; 
      }
    } 
    //Check for Nan
    if (checkNan(u) || checkNan(v) || checkNan(rho)) {
      std::cout << "Iteration : " << t+1 << " ...exiting\n";
      return 1;
    }
    if (!((t+1)%printFrequency)) {
      std::cout << "Iteration : " << t+1 << "\t";
      int index1 = nx-1;
      int index2 = ny/2-1;
      std::cout << rho[index1][index2]<< "\t" << u[index1][index2] << "\t"
                << v[index1][index2] << "\n";
    }
  }
  dumpTecplot(u, v);
  dumpTecplotMag(u, v);
  dumpBinary(u, v);
  return 0;
}
void dumpTecplot(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1.open("fpcTecPlot.dat");
  f1 << std::scientific;
  f1.precision(4);
  f1 << "TITLE = Re=" << Re << "100 Flow Past Cylinder" << std::endl \
     << "VARIABLES = X, Y, U, V" << std::endl;
  f1 << "Zone I =" << nx << " J=" << ny << "F=POINT" << std::endl;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int position = std::pow(i-circleX, 2)+std::pow(j-circleY, 2);
      if ( position <= circleR2) {
        double v = 0.0;
        f1 << i << "\t" << j << "\t" << v << "\t" << v << "\n";
      } else {
        f1 << i << "\t" << j << "\t" << u[i][j] << "\t" << v[i][j] << "\n";
      }
    }
  }
  f1.close();
}
void dumpTecplotMag(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1.open("fpcMagTecPlot.dat");
  f1 << std::scientific;
  f1.precision(4);
  f1 << "TITLE = Re=" << Re << "100 Flow Past Cylinder" << std::endl \
     << "VARIABLES = X, Y, U" << std::endl;
  f1 << "Zone I =" << nx << " J=" << ny << "F=POINT" << std::endl;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int position = std::pow(i-circleX, 2)+std::pow(j-circleY, 2);
      if ( position <= circleR2) {
        double v = 0.0;
        f1 << i << "\t" << j << "\t" << v << "\n";
      } else {
        f1 << i << "\t" << j << "\t" << std::sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]) << "\n";
      }
    }
  }
  f1.close();
}
void dumpBinary(double u[][ny], double v[][ny]) {
  std::ofstream f1;
  f1.open("fpcBin.dat", std::ios::binary);
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double x = (i+0.5)/nx;
      double y = (j+0.5)/ny;
      double vmag = std::sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
      int position = std::pow(i-circleX, 2)+std::pow(j-circleY, 2);
      if ( position <= circleR2) {
        vmag = 0.0;
      } 
      f1.write((char*)&x, sizeof(double));
      f1.write((char*)&y, sizeof(double));
      f1.write((char*)&vmag, sizeof(double));
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
