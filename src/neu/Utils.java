package neu;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;



public class Utils {
 
  public static float[] pyrimidalCellFiring(
    Sampling st, float ij, int ns) 
  {
    //float m0 = 0.001f;
    //float h0 = 0.999f;
    //float n0 = 0.001f;
    //float s0 = 0.009f;
    //float c0 = 0.007f;
    //float q0 = 0.010f;
    //float vd0 = v0-4.5f;
    //float vs0 = v0-4.6f;
    //if (ca0==null) float ca0 = 0.2f;
    // ion channel kinetics initial conditions
    float n0 = 0.0003f;
    float m0 = 0.0011f;
    float h0 = 0.9998f;
    float v0 = -65.0f;
    int nt = st.getCount();
    int isi = nt/ns;
    float[] v = new float[nt];
    double dt = st.getDelta();
    Sampling stt = new Sampling(isi,dt,st.getFirst());
    for (int i=0; i<ns; ++i) { 
      float[] vt = rk4HH(stt,v0,n0,m0,h0,ij);
      copy(isi,0,vt,i*isi,v);
      ++i;
      if (i>=ns) break;
      stt = new Sampling(isi,dt,i*isi*dt);
      float vl = vt[isi-1];
      float[] vi = rk4L(stt,vl,ij);
      copy(isi,0,vi,i*isi,v);
      stt = new Sampling(isi,dt,i*isi*dt);
    }
    //float[] v = rk4HH(st,v0,n0,m0,h0,ij);
    return v;
  }

  public static float[] rk4HH(
    Sampling st, float Vm, float n0, 
    float m0, float h0, float ij) 
  {
    int nt = st.getCount();
    int ntm = nt-1;
    float dt = (float)st.getDelta();
    float[] f = new float[nt];
    f[0] = Vm;
    FinalParams p = new FinalParams();
    for (int i=0; i<ntm; ++i) {
      float fi = f[i];
      float mfi = -fi;
      float t = (float)st.getValue(i);
      //f[i+1] = f[i] + dt*vprime(f[i],m0,n0,h0,inj); //euler
      float v1 = dt*p.vprime(t,fi,m0,n0,h0,ij,p.Ica(fi));
      float v2 = dt*p.vprime(t+dt/2.0f,fi+v1/2.0f,m0,n0,h0,ij,p.Ica(fi));
      float v3 = dt*p.vprime(t+dt/2.0f,fi+v2/2.0f,m0,n0,h0,ij,p.Ica(fi));
      float v4 = dt*p.vprime(t+dt,fi+v3,m0,n0,h0,ij,p.Ica(fi));
      f[i+1] = f[i] + (v1+2*v2+2*v3+v4)/6.0f;
      v1 = dt*p.nprime(t,n0,mfi); 
      v2 = dt*p.nprime(t+dt/2.0f,n0,mfi+v1/2.0f); 
      v3 = dt*p.nprime(t+dt/2.0f,n0,mfi+v2/2.0f); 
      v4 = dt*p.nprime(t+dt,n0,mfi); 
      n0 += (v1+2*v2+2*v3+v4)/6.0f;
      v1 = dt*p.mprime(t,m0,mfi); 
      v2 = dt*p.mprime(t+dt/2.0f,m0,mfi+v1/2.0f); 
      v3 = dt*p.mprime(t+dt/2.0f,m0,mfi+v2/2.0f); 
      v4 = dt*p.mprime(t+dt,m0,mfi); 
      m0 += (v1+2*v2+2*v3+v4)/6.0f;
      v1 = dt*p.hprime(t,h0,mfi); 
      v2 = dt*p.hprime(t+dt/2.0f,h0,mfi+v1/2.0f); 
      v3 = dt*p.hprime(t+dt/2.0f,h0,mfi+v2/2.0f); 
      v4 = dt*p.hprime(t+dt,h0,mfi); 
      h0 += (v1+2*v2+2*v3+v4)/6.0f;
    }
    return f;
  }

  public static float[] rk4L(Sampling st, float v0, float ij) {
    int nt = st.getCount();
    int ntm = nt-1;
    float dt = (float)st.getDelta();
    float[] f = new float[nt];
    f[0] = v0;
    FinalParams p = new FinalParams();
    for (int i=0; i<ntm; ++i) {
      float t = (float)st.getValue(i);
      float v1 = dt*p.lprime(t,f[i],ij);
      float v2 = dt*p.lprime(t+dt/2.0f,f[i]+v1/2.0f,ij);
      float v3 = dt*p.lprime(t+dt/2.0f,f[i]+v2/2.0f,ij);
      float v4 = dt*p.lprime(t+dt,f[i]+v3,ij);
      f[i+1] = f[i] + (v1+2*v2+2*v3+v4)/6.0f;
    }
    return f;
  }

  public static float[] rk4dca(Sampling st, float c0, float tau, float[] ica) {
    int nt = st.getCount();
    int ntm = nt-1;
    float dt = (float)st.getDelta();
    float[] c = new float[nt];
    c[0] = c0;
    float l = 1.0f;
    for (int i=0; i<ntm; ++i) {
      float t = (float)st.getValue(i);
      float v1 = dt*dcdt(t,c[i],ica[i],tau); 
      float v2 = dt*dcdt(t+dt/2.0f,c[i]+v1/2.0f,ica[i],tau);
      float v3 = dt*dcdt(t+dt/2.0f,c[i]+v2/2.0f,ica[i],tau);
      float v4 = dt*dcdt(t+dt,c[i]+v3,ica[i],tau);
      c[i+1] = c[i] + (v1+2*v2+2*v3+v4)/6.0f;
    }
    return c;
  }
  public static float dcdt(float t, float ca, float ica, float tau) {
    return ica-ca/tau;
  }


  public static float[] rk4dw(
    Sampling st, float[] ca, float[] eta, 
    float[] omega, float w0) 
  {
    int nt = st.getCount();
    int ntm = nt-1;
    float dt = (float)st.getDelta();
    float[] w = new float[nt];
    w[0] = w0;
    float l = 1.0f;
    for (int i=0; i<ntm; ++i) {
      float t = (float)st.getValue(i);
      float v1 = dt*dwdt(t,w[i],eta[i],omega[i],l); 
      float v2 = dt*dwdt(t+dt/2.0f,w[i]+v1/2.0f,eta[i],omega[i],l);
      float v3 = dt*dwdt(t+dt/2.0f,w[i]+v2/2.0f,eta[i],omega[i],l);
      float v4 = dt*dwdt(t+dt,w[i]+v3,eta[i],omega[i],l);
      w[i+1] = w[i] + (v1+2*v2+2*v3+v4)/6.0f;
    }
    return w;
  }
  /**
   * Calcium control hypothesis ODE.
   * @param t     time
   * @param eta   Ca dependent learning rate (constant)
   * @param omega Ca dependent sign magnitude (constant)
   * @param w     synaptic plasticity weight 
   * @param l     decay constant 
   */
  public static float dwdt(float t, float w, float eta, float omega, float l) {
    return eta*(omega-w*l);
  }


  private static class FinalParams {
    double C   = 1.0;
    double gk  = 0.360*100;//36.0;
    double gna = 2.000*100;//200.0;
    double gl  = 0.3;
    double gca = 0.100*100;
    double gld = 0.030*100;
    double Ek  =-80.0;
    double Ena = 50.0;
    double El  =-65.0;
    double Eca = 75.0;
    public float mprime(float t, float m, float V) {
      double am = 0.1*(V+25.0)/(exp((V+25.0)/10.0)-1.0);
      double bm = 4.0*exp(V/18.0);
      return (float)((1.0-m)*am-m*bm);
    }
    public float nprime(float t, float n, float V) {
      double an = 0.01*(V+10.0)/(exp((V+10.0)/10.0)-1.0);
      double bn = 0.125*exp(V/80.0);
      return (float)((1.0-n)*an-n*bn);
    }
    public float hprime(float t, float h, float V) {
      double ah = 0.07*exp(V/20.0);
      double bh = 1.0/(exp((V+30.0)/10.0)+1.0);
      return (float)((1.0-h)*ah-h*bh);
    }
    public float lprime(float t, float V, float ij) {
      return (float)((1.0/C)*(ij-gld*(V-El)));
    }
    public float Ica(float v) {
      double mgo = 1.0;
      double vr = 30.77;
      double bv = 1.0/(1.0+mgo*exp(-0.062*v)/3.57);
      return (float)(bv*(v-vr));
    }
    public float vprime(
      float t, float Vm, float m, 
      float n, float h, float ij, float ica) 
    {
      double dV = (-1.0/C)*(gk*(n*n*n*n)*(Vm-Ek)
                           +gna*(m*m*m)*h*(Vm-Ena)
                           +gl*(Vm-El)-ij);
      return (float)dV;
    }
  }


//  private static class Traub {
//    // maximal conductances in mS/cm^2
//    float gl = 0.1f; 
//    float gna = 30.0f; 
//    float gca = 10.0f;
//    float gkdr = 15.0f;
//    float gkahp = 0.8f;
//    float gkc = 15.0f;
//    float gnmda = 0.0f;
//    float gampa = 0.0f;
//    // reverse potentials in mV
//    float ena = 120.0f;
//    float eca = 140.0f;
//    float ek  = -15.0f;
//    float el  = 0.0f;
//    float esy = 60.0f;
//    // applied currents in uA/cm^2
//    float Is = -0.5f;
//    float Id = 0.0f;
//    // coupling parameters
//    float gc = 2.1; // mS/cm^2
//    float p  = 0.5;
//    float C = 3; // uF/cm^2
//    public static float chi(float ca) {
//      return (float)min(ca/250.0,1.0);
//    }
//    public static float mprime(float t, float m, float vs) {
//      double am = 0.32*(13.1-vs)/(exp((13.1-vs)/4.0)-1.0);
//      double bm = 0.28*(vs-40.1)/(exp((vs-40.1)/5.0)-1.0);
//      return (float)(1.0-m)*am-m*bm;
//    }
//    public static float nprime(float t, float n, float vs) {
//      double an = 0.016*(35.1-vs)/(exp((35.1-vs)/5.0)-1.0);
//      double bn = 0.25*exp(0.5-0.25*vs);
//      return (float)(1.0-n)*an-n*bn;
//    }
//    public static float hprime(float t, float h, float vs) {
//      double ah = 0.128*exp((17.0-vs)/18.0);
//      double bh = 4.0/(1.0+exp(40.0-vs)/5.0);
//      return (float)(1.0-h)*ah-h*bh;
//    }
//    public static float sprime(float t, float s, float vd) {
//      double as = 1.6/(1.0+exp(-0.072*(vd-65.0)));
//      double bs = 0.02*(vd-51.1)/(exp((vd-51.1)/5.0)-1.0);
//      return (float)(1.0-s)*as-s*bs;
//    }
//    public static float cprime(float t, float c, float vd) {
//      // for vd \leq 50
//      if (vd<=50.0) {
//        double ac = (exp((vd-10.0)/11.0)-exp(vd-6.5)/27.0)/18.975;
//        double bc = 2.0*exp((6.5-vd)/27.0)-ac;
//      }
//      if (vd>50.0) {
//      // for vd \gt 50
//        double ac = 2.0*exp((6.5-vd)/27.0);
//        double bc = 0.0;
//      }
//      return (float)(1.0-c)*ac-c*bc;
//    }
//    public static float qprime(float t, float q, float ca) {
//      double aq = min((0.00002)*ca,0.01);
//      double bq = 0.001;
//      return (float)(1.0-q)*aq-q*bq;
//    }
//  }
//   private static class Warman {
//    // maximal conductances in mS/cm^2
//    float gl = 0.1f; 
//    float gna = 150.0f; 
//    float gca = 10.0f;
//    float gkdr = 15.0f;
//    float gkahp = 0.8f;
//    float gkc = 15.0f;
//    float gnmda = 0.0f;
//    float gampa = 0.0f;
//    // reverse potentials in mV
//    float ena = 65.0f;
//    float eca = -25.8*log(1.0/2000.0);
//    float ek  = -80.0f;
//    float el  = 0.0f;
//    float esy = 60.0f;
//    // applied currents in uA/cm^2
//    float Is = -0.5f;
//    float Id = 0.0f;
//    // coupling parameters
//    float gc = 2.1; // mS/cm^2
//    float p  = 0.5;
//    float C = 3; // uF/cm^2
//    public static float chi(float ca) {
//      return (float)min(ca/250.0,1.0);
//    }
//    public static float mprime(float t, float m, float v) {
//      double am = -1.74*(v-11.0)/(exp((v-11.0)/-12.94)-1.0);
//      double bm = 0.06*(v-5.9)/(exp((v-5.9)/4.47)-1.0);
//      return (float)(1.0-m)*am-m*bm;
//    }
//    public static float hprime(float t, float h, float v) {
//      double ah = 3.0/exp((v+80.0)/10.0);
//      double bh = 12.0/(exp((v-77.0)/-27.0)+1.0);
//      return (float)(1.0-h)*ah-h*bh;
//    }
//    public static float nprime(float t, float n, float v) {
//      double an = 0.016*(35.1-vs)/(exp((35.1-vs)/5.0)-1.0);
//      double bn = 0.25*exp(0.5-0.25*vs);
//      return (float)(1.0-n)*an-n*bn;
//    }
//    public static float sprime(float t, float s, float v) {
//      double as = -0.16*(v+26.0)/(exp((v+26.0)/-4.5)-1.0);
//      double bs =  0.04*(v+12.0)/(exp((v+12.0)/10.0)-1.0);
//      return (float)(1.0-s)*as-s*bs;
//    }
//    public static float rprime(float t, float r, float v) {
//      double ar = 2.0/(exp((v+94.0)/10.0));
//      double br = 8.0/(exp((v-68.0)/-27.0)+1.0);
//      return (float)(1.0-r)*ar-r*br;
//    }
//    public static float qprime(float t, float q, float ca) {
//      double aq = min((0.00002)*ca,0.01);
//      double bq = 0.001;
//      return (float)(1.0-q)*aq-q*bq;
//    }
//}



}
