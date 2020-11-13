//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk.cpp
//  \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//  spherical polar coordinates.  Initial conditions are in vertical hydrostatic eqm.
//  Editted for the purpose of studying the vortex formation

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
void VelProfileCyl(const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3);
void Planet_and_Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// declaration of hitory output file
Real hist_rp(MeshBlock *pmb, int iout);
Real hist_phip(MeshBlock *pmb, int iout);
Real hist_acc1p(MeshBlock *pmb, int iout);
Real hist_acc2p(MeshBlock *pmb, int iout);
Real hist_massp(MeshBlock *pmb, int iout);
Real hist_rx(MeshBlock *pmb, int iout);
Real hist_ry(MeshBlock *pmb, int iout);
Real hist_vx(MeshBlock *pmb, int iout);
Real hist_vy(MeshBlock *pmb, int iout);
Real hist_oldrx(MeshBlock *pmb, int iout);
Real hist_oldry(MeshBlock *pmb, int iout);
Real hist_oldvx(MeshBlock *pmb, int iout);
Real hist_oldvy(MeshBlock *pmb, int iout);
Real hist_ax(MeshBlock *pmb, int iout);
Real hist_ay(MeshBlock *pmb, int iout);
Real hist_dt(MeshBlock *pmb, int iout);

// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas, mass_S;
Real r_P, phi_P, vr_P, vphi_P, final_mass_P, tau;
Real soft, soft_in;
Real acc1_P, acc2_P;
Real old_acc1_P, old_acc2_P;
Real acc1, acc2, src1, src2;
Real mass_P;

Real r_x, r_y;
Real v_x, v_y;
Real old_r_x, old_r_y;
Real old_v_x, old_v_y;
Real old_a_x, old_a_y;
Real a_x, a_y;

int dt_step, rk, rk_s;
//Real omega_P;

} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0); // mass of Sun, G=1 in here
  r0 = pin->GetOrAddReal("problem","r0",1.0); // r0 is the unit of r(AU)
  mass_S = 1.0;

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem", "dslope", 1.5);

  // Get parameters of initial pressure and cooling parameters
  // Note that non barotropic is set as 1 in default
  gamma_gas = pin->GetReal("hydro","gamma");
  // p0_over_r0 = initial pressure/density without over gamma(using ideal gas law), = SQR(sound speed)
  p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025); // SQR(sound speed)
  pslope = pin->GetOrAddReal("problem", "pslope", 2.5);

  // Planet initial parameter
  r_P = 1.0; //initial position of r in AU
  phi_P = 0.0; //initial position of phi
  final_mass_P = pin->GetOrAddReal("problem", "final_mass_P", 1.0)/1047.35; //in the unit of solar mass(from jupiter mass)
  vr_P = 0.0;
  vphi_P = sqrt(gm0/r_P);
  acc1_P = -gm0/pow(r_P,2);
  acc2_P = 0;

  dt_step = 1;
  rk = pin->GetOrAddReal("problem", "rk", 1);
  rk_s = pin->GetOrAddReal("problem", "rk_s", 1);

  //omega_P = 1;

//  old_r_x = r_P;
//  old_r_y = 0;
//  r_x = old_r_x;
//  r_y = old_r_y;
//
//  old_v_x = 0;
//  old_v_y = 1;
//  v_x = old_v_x;
//  v_y = old_v_y;
//
//  a_x = acc1_P;
//  a_y = 0;

  tau = pin->GetOrAddReal("problem", "tau", 0.01);

  // Enroll the source function
  EnrollUserExplicitSourceFunction(Planet_and_Cooling);

  // user defined output in histry file
  AllocateUserHistoryOutput(9);
  EnrollUserHistoryOutput(0, hist_rp, "r_P");
  EnrollUserHistoryOutput(1, hist_phip, "phi_P");
  EnrollUserHistoryOutput(2, hist_acc1p, "acc1p");
  EnrollUserHistoryOutput(3, hist_acc2p, "acc2p");
  EnrollUserHistoryOutput(4, hist_massp, "mass_P");
  EnrollUserHistoryOutput(5, hist_rx, "r_x");
  EnrollUserHistoryOutput(6, hist_ry, "r_y");
  EnrollUserHistoryOutput(7, hist_vx, "v_x");
  EnrollUserHistoryOutput(8, hist_vy, "v_y");
//  EnrollUserHistoryOutput(9, hist_oldrx, "old_rx");
//  EnrollUserHistoryOutput(10, hist_oldry, "old_ry");
//  EnrollUserHistoryOutput(11, hist_oldvx, "old_vx");
//  EnrollUserHistoryOutput(12, hist_oldvy, "old_vy");
//  EnrollUserHistoryOutput(13, hist_dt, "h_dt");
//  EnrollUserHistoryOutput(14, hist_ax, "ax");
//  EnrollUserHistoryOutput(15, hist_ay, "ay");

  
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);

        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;

        Real p_over_r = PoverR(rad,phi,z);
        phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }
  }

  return;
}

//


//void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//{
//    AllocateUserOutputVariables(6);
//    return;
//}
//
//void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
//    for(int k=ks; k<=ke; k++) {
//        for(int j=js; j<=je; j++) {
//            for(int i=is; i<=ie; i++) {
//            user_out_var(0,k,j,i) = r_P;
//            user_out_var(1,k,j,i) = phi_P;
//            user_out_var(2,k,j,i) = acc1_P;
//            user_out_var(3,k,j,i) = acc2_P;
//            user_out_var(4,k,j,i) = src1;
//            user_out_var(5,k,j,i) = src2;
//          }
//        }
//      }
//}


namespace {
//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  den = rho0*std::pow(rad/r0 , -dslope);
  return den;
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, -1)/gamma_gas;
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

void VelProfileCyl(const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3){
  Real vel;
  vel = std::sqrt(gm0/rad + (p0_over_r0/gamma_gas)*(-pslope)*std::pow(rad/r0 , -1));
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    v1=0.0;
    v2=vel;
    v3=0.0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    v1=0.0;
    v2=0.0;
    v3=vel;
  }
  return;
}

void Planet_and_Cooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  // MeshBlock *pmb = pmy_hydro_->pmy_block;
  // pmb is already assigned before calling this source term function
  Real src1_P, src2_P;
  Real accx, accy, accx_P, accy_P;

  Real r_P_cube = std::pow(1.0,3); //initial r_P in cube
  Real time_P = 2*M_PI*std::pow(r_P_cube/gm0 , 0.5); //period of planet
  AthenaArray<Real> vol(pmb->ncells1);

  if (time < time_P) {
      mass_P = final_mass_P*std::pow(std::sin((PI/2)*(time/time_P)),2);
  }
  else {
      mass_P = final_mass_P;
  }

  acc1_P = -gm0/pow(r_P,2);
  acc2_P = 0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
// #pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
// #pragma simd
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real rad, phi, z;
        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k);
        Real theta = phi - phi_P; // rel. location of the planet
        //calculate the force on disk element
        soft = 0.6*r_P*std::pow(mass_P/(3*mass_S),1.0/3); // soft = k*R_H
        //dist is the distance between planet and each disk element
        Real dist = std::sqrt(rad*rad + r_P*r_P - 2.0*rad*r_P*std::cos(theta));
        // distance to the planet with soften parameter in the G-potential
        Real dist_soft = std::sqrt(rad*rad + r_P*r_P - 2.0*rad*r_P*std::cos(theta) + std::pow(soft,2));
        // angle opposite to rad(radius of disk element w.r.p. to sun)
        Real alpha = std::acos((r_P*r_P+dist*dist-rad*rad)/(2.0*dist*r_P));
        // accound for the direction below x-axis for accy
        if (theta > PI){
            alpha = -alpha;
        }
        // acceleration in r and phi
        // projection on the vector r
        // for r so far away, force at x-dir is +ve, for r very small, force at x-dir is-ve
        // vector_dist dot unit_vector x
        accx = mass_P/std::pow(dist_soft,3.0)*(dist*cos(2*M_PI-alpha));
        // vector_dist dot unit_vector y
        accy = mass_P/std::pow(dist_soft,3.0)*(dist*cos(1.5*M_PI-alpha));
        // rotate the force back to the cylindrical coordinate
        acc1 = cos(theta)*accx + sin(theta)*accy;
        acc2 = -sin(theta)*accx + cos(theta)*accy;

        // sources to the radial and angular momentum eqs
        src1 = dt*den*acc1;
        src2 = dt*den*acc2;
        // update
        cons(IM1,k,j,i) += src1;
        cons(IM2,k,j,i) += src2;
        cons(IEN,k,j,i) += 0.5*(SQR(src1)+SQR(src2))/den;

        // calculate the force on planet
        //notice that the direction of vector dist is reversed, so we don't need a -ve sign in front of acc_P
        soft_in = 1 - std::exp(-dist*dist/std::pow(soft,2));
        accx_P = -den*vol(i)/std::pow(dist_soft,3.0)*(dist*cos(2*M_PI-alpha));
        accy_P = -den*vol(i)/std::pow(dist_soft,3.0)*(dist*cos(1.5*M_PI-alpha)); //*vol(i)
        acc1_P += soft_in*accx_P;
        acc2_P += soft_in*accy_P;
        // short description: the above calculate the acceleration by the planet on the elements of the disk
        // and add the acceration according to the direction of vector r and vector phi of each element
        // note that the calculation is using theta, which means we use the planet as a reference

      }
    }
  }

  //To Dr. Leung, if you want to test the if statement
  //please unhide this part and  comment-out the LF integrator
//  phi_P = omega_P*time;
//  if (phi_P > TWO_PI){
//      phi_P = phi_P - TWO_PI;
//  }

  //LF
  //For the first interaction, the above is computing a_i
  //in order to compute v_i+1, we require a_i+1, which is computed in the next iteration
  //that's why we save the i-iteraction's acc_P as old_acc_P for the next iteration's v_i+1

  //LF integrator
  // If Rk1, dt_step % 1, rk_s = 0
  // If RK2, dt_step % 2, rk_s = 1
  // If RK3, dt_step % 3, rk_s = 1
  if (dt_step % rk == rk_s) {
      old_r_x = r_P*cos(phi_P);
      old_r_y = r_P*sin(phi_P);
      old_v_x = vr_P*cos(phi_P) - vphi_P*sin(phi_P);
      old_v_y = vr_P*sin(phi_P) + vphi_P*cos(phi_P);

      a_x = acc1_P*cos(phi_P) - acc2_P*sin(phi_P);
      a_y = acc1_P*sin(phi_P) + acc2_P*cos(phi_P);

      // This is the calculation of v_i, since we need a_(i+1) for v_i, so we hold the calculation of v_i
      // at the initial time.
      // x_(i+1) = x_i +v_i*dt + 0.5*a_i*dt*dt
      // v_(i+1) = v_i + 0.5*(a_i+a_(i+1))*dt*dt

      // notice that for x1 we need v_0 only
      // so we don't need to compute v at this moment
      if (time > dt){
          v_x = old_v_x + 0.5*(old_a_x +  a_x)*(dt);
          v_y = old_v_y + 0.5*(old_a_y +  a_y)*(dt);
      } else {
          v_x = old_v_x;
          v_y = old_v_y;
      }

      r_x = old_r_x + v_x*(dt) + 0.5*a_x*(dt)*(dt);
      r_y = old_r_y + v_y*(dt) + 0.5*a_y*(dt)*(dt);

      old_a_x = a_x;
      old_a_y = a_y;

      r_P = sqrt(r_x*r_x + r_y*r_y);
      phi_P = std::atan2(r_y, r_x);

      if (phi_P < 0) {
          phi_P += 2*M_PI;
      }

      vr_P = v_x*cos(phi_P) + v_y*sin(phi_P);
      vphi_P = -v_x*sin(phi_P) + v_y*cos(phi_P);
  }

  dt_step += 1;

  // Rmk: the above leapfrog integrator is using a fixed coordinate setted by athena++,
  // that's why we rotate the system using phi_P

  // Euler integrator
//  old_r_x = r_P*cos(phi_P);
//  old_r_y = r_P*sin(phi_P);
//
//  old_v_x = vr_P*cos(phi_P) - vphi_P*sin(phi_P);
//  old_v_y = vr_P*sin(phi_P) + vphi_P*cos(phi_P);
//
//  a_x = acc1_P*cos(phi_P) - acc2_P*sin(phi_P);
//  a_y = acc1_P*sin(phi_P) + acc2_P*cos(phi_P);
//
//  v_x = old_v_x + a_x*(dt*2/3);
//  v_y = old_v_y + a_y*(dt*2/3);
//
//  r_x = old_r_x + old_v_x*(dt*2/3);
//  r_y = old_r_y + old_v_y*(dt*2/3);
//
//  r_P = sqrt(r_x*r_x + r_y*r_y);
//  phi_P = std::atan2(r_y, r_x);
//
//  vr_P = v_x*cos(phi_P) + v_y*sin(phi_P);
//  vphi_P = -v_x*sin(phi_P) + v_y*cos(phi_P);

  // Euler-Cromer integrator
//  old_r_x = r_P*cos(phi_P);
//  old_r_y = r_P*sin(phi_P);
//
//  old_v_x = vr_P*cos(phi_P) - vphi_P*sin(phi_P);
//  old_v_y = vr_P*sin(phi_P) + vphi_P*cos(phi_P);
//
//  a_x = acc1_P*cos(phi_P) - acc2_P*sin(phi_P);
//  a_y = acc1_P*sin(phi_P) + acc2_P*cos(phi_P);
//
//  v_x = old_v_x + a_x*(dt*2/3);
//  v_y = old_v_y + a_y*(dt*2/3);
//
//  r_x = old_r_x + v_x*(dt*2/3);
//  r_y = old_r_y + v_y*(dt*2/3);
//
//  r_P = sqrt(r_x*r_x + r_y*r_y);
//  phi_P = std::atan2(r_y, r_x);
//
//  vr_P = v_x*cos(phi_P) + v_y*sin(phi_P);
//  vphi_P = -v_x*sin(phi_P) + v_y*cos(phi_P);

  //since we CANNOT enroll two source term function, so we combine the source functions together
  //cooling function

  Real omega;
  Real tau_r;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
          pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
          for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real rad, phi, z;
              GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k);
              Real den = prim(IDN,k,j,i);
              Real poverr;
              omega = prim(IVY,k,j,i)/rad;
              tau_r = 2*M_PI*tau/omega;
              poverr = p0_over_r0*std::pow(rad/r0, -1)/gamma_gas;
              Real new_pres = prim(IPR,k,j,i) - dt/tau_r*(prim(IPR,k,j,i) - poverr*den);
              cons(IEN,k,j,i) += (new_pres - prim(IPR,k,j,i))/(gamma_gas-1);
              // cons(IEN,k,j,i) = new_pres/(gamma_gas-1);
              // cons(IEN,k,j,i) += 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i)))/den;

              // cons(IEN,k,j,i) += 0.5*den*(SQR(prim(IVX,k,j,i))+SQR(prim(IVY,k,j,i)));

              // cons(IEN,k,j,i) += (-dt/tau_r*(cons(IEN,k,j,i)*(gamma_gas-1)-poverr*den))/(gamma_gas-1);
          }
      }
  }

  return;
}

// hitory output file function
Real hist_rp(MeshBlock *pmb, int iout){
  return r_P;
}

Real hist_phip(MeshBlock *pmb, int iout){
  return phi_P;
}

Real hist_acc1p(MeshBlock *pmb, int iout){
  return acc1_P;
}

Real hist_acc2p(MeshBlock *pmb, int iout){
  return acc2_P;
}

Real hist_massp(MeshBlock *pmb, int iout){
  return mass_P;
}

Real hist_rx(MeshBlock *pmb, int iout){
  return r_x;
}

Real hist_ry(MeshBlock *pmb, int iout){
  return r_y;
}

Real hist_vx(MeshBlock *pmb, int iout){
  return v_x;
}

Real hist_vy(MeshBlock *pmb, int iout){
  return v_y;
}

Real hist_oldrx(MeshBlock *pmb, int iout){
  return old_r_x;
}

Real hist_oldry(MeshBlock *pmb, int iout){
  return old_r_y;
}

Real hist_oldvx(MeshBlock *pmb, int iout){
  return old_v_x;
}

Real hist_oldvy(MeshBlock *pmb, int iout){
  return old_v_y;
}

Real hist_ax(MeshBlock *pmb, int iout){
  return a_x;
}

Real hist_ay(MeshBlock *pmb, int iout){
  return a_y;
}

} // namespace


