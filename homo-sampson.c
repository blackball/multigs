/**
 * Caculate Sampson error for homography estimation.
 *
 * The codes were adapted from Homest 1.4
 */

#include <math.h>
static void
homo_sampson(const double m1[2], const double m2[2],const double h[9],double err[1]) {
        double t1;
        double t10;
        double t100;
        double t104;
        double t108;
        double t112;
        double t118;
        double t12;
        double t122;
        double t125;
        double t126;
        double t129;
        double t13;
        double t139;
        double t14;
        double t141;
        double t144;
        double t15;
        double t150;
        double t153;
        double t161;
        double t167;
        double t17;
        double t174;
        double t18;
        double t19;
        double t193;
        double t199;
        double t2;
        double t20;
        double t201;
        double t202;
        double t21;
        double t213;
        double t219;
        double t22;
        double t220;
        double t222;
        double t225;
        double t23;
        double t236;
        double t24;
        double t243;
        double t250;
        double t253;
        double t26;
        double t260;
        double t27;
        double t271;
        double t273;
        double t28;
        double t29;
        double t296;
        double t3;
        double t30;
        double t303;
        double t31;
        double t317;
        double t33;
        double t331;
        double t335;
        double t339;
        double t34;
        double t342;
        double t345;
        double t35;
        double t350;
        double t354;
        double t36;
        double t361;
        double t365;
        double t37;
        double t374;
        double t39;
        double t4;
        double t40;
        double t41;
        double t42;
        double t43;
        double t44;
        double t45;
        double t46;
        double t47;
        double t49;
        double t51;
        double t57;
        double t6;
        double t65;
        double t66;
        double t68;
        double t69;
        double t7;
        double t72;
        double t78;
        double t8;
        double t86;
        double t87;
        double t90;
        double t95;
        {
                t1 = m2[0];
                t2 = h[6];
                t3 = t2*t1;
                t4 = m1[0];
                t6 = h[7];
                t7 = t1*t6;
                t8 = m1[1];
                t10 = h[8];
                t12 = h[0];
                t13 = t12*t4;
                t14 = h[1];
                t15 = t14*t8;
                t17 = t3*t4+t7*t8+t1*t10-t13-t15-h[2];
                t18 = m2[1];
                t19 = t18*t18;
                t20 = t2*t2;
                t21 = t19*t20;
                t22 = t18*t2;
                t23 = h[3];
                t24 = t23*t22;
                t26 = t23*t23;
                t27 = t6*t6;
                t28 = t19*t27;
                t29 = t18*t6;
                t30 = h[4];
                t31 = t29*t30;
                t33 = t30*t30;
                t34 = t4*t4;
                t35 = t20*t34;
                t36 = t2*t4;
                t37 = t6*t8;
                t39 = 2.0*t36*t37;
                t40 = t36*t10;
                t41 = 2.0*t40;
                t42 = t8*t8;
                t43 = t42*t27;
                t44 = t37*t10;
                t45 = 2.0*t44;
                t46 = t10*t10;
                t47 = t21-2.0*t24+t26+t28-2.0*t31+t33+t35+t39+t41+t43+t45+t46;
                t49 = t12*t12;
                t51 = t6*t30;
                t57 = t20*t2;
                t65 = t1*t1;
                t66 = t65*t20;
                t68 = t65*t57;
                t69 = t4*t10;
                t72 = t2*t49;
                t78 = t27*t6;
                t86 = t65*t78;
                t87 = t8*t10;
                t90 = t65*t27;
                t95 = -2.0*t49*t18*t51-2.0*t3*t12*t46-2.0*t1*t57*t12*t34-2.0*t3*t12*t33+t66
                        *t43+2.0*t68*t69+2.0*t72*t69-2.0*t7*t14*t46-2.0*t1*t78*t14*t42-2.0*t7*t14*t26+
                        2.0*t86*t87+t90*t35+2.0*t49*t6*t87;
                t100 = t14*t14;
                t104 = t100*t2;
                t108 = t2*t23;
                t112 = t78*t42*t8;
                t118 = t57*t34*t4;
                t122 = t10*t26;
                t125 = t57*t4;
                t126 = t10*t19;
                t129 = t78*t8;
                t139 = -2.0*t57*t34*t18*t23+2.0*t100*t6*t87+2.0*t104*t69-2.0*t100*t18*t108+
                        4.0*t36*t112+6.0*t43*t35+4.0*t118*t37+t35*t28+2.0*t36*t122+2.0*t125*t126+2.0*
                        t129*t126+2.0*t37*t122-2.0*t78*t42*t18*t30+t43*t21;
                t141 = t10*t33;
                t144 = t46*t18;
                t150 = t46*t19;
                t153 = t46*t10;
                t161 = t27*t27;
                t167 = 2.0*t36*t141-2.0*t144*t108+2.0*t37*t141+t66*t33+t150*t27+t150*t20+
                        4.0*t37*t153+6.0*t43*t46+4.0*t112*t10+t43*t33+t161*t42*t19+t43*t26+4.0*t36*t153
                        ;
                t174 = t20*t20;
                t193 = 6.0*t35*t46+4.0*t10*t118+t35*t33+t35*t26+t174*t34*t19+t100*t27*t42+
                        t100*t20*t34+t100*t19*t20+t90*t46+t65*t161*t42+t90*t26+t49*t27*t42+t49*t20*t34+
                        t49*t19*t27;
                t199 = t34*t34;
                t201 = t12*t23;
                t202 = t14*t30;
                t213 = t42*t42;
                t219 = t66*t46+t100*t26+t46*t100+t174*t199-2.0*t201*t202-2.0*t144*t51+t46*
                        t26+t65*t174*t34+t49*t33+t49*t46+t46*t33+t161*t213-2.0*t7*t14*t20*t34;
                t220 = t1*t27;
                t222 = t36*t8;
                t225 = t7*t14;
                t236 = t4*t6*t8;
                t243 = t3*t12;
                t250 = t46*t46;
                t253 = t1*t20;
                t260 = -4.0*t220*t14*t222-4.0*t225*t40-4.0*t220*t15*t10+2.0*t90*t40+2.0*
                        t225*t24+2.0*t72*t236-2.0*t3*t12*t27*t42-4.0*t243*t44+2.0*t66*t44+2.0*t243*t31+
                        t250+2.0*t68*t236-4.0*t253*t12*t236-4.0*t253*t13*t10;
                t271 = t4*t20;
                t273 = t8*t18;
                t296 = t10*t18;
                t303 = 2.0*t104*t236-2.0*t35*t31+12.0*t35*t44+2.0*t125*t37*t19-4.0*t271*t6*
                        t273*t23+2.0*t36*t37*t26+2.0*t36*t129*t19-4.0*t36*t27*t273*t30+2.0*t36*t37*t33+
                        12.0*t36*t43*t10+12.0*t36*t37*t46-4.0*t271*t296*t23+2.0*t36*t126*t27;
                t317 = t18*t14;
                t331 = t14*t2;
                t335 = t12*t18;
                t339 = t220*t18;
                t342 = t7*t30;
                t345 = t317*t6;
                t350 = -4.0*t31*t40-2.0*t43*t24+2.0*t37*t126*t20-4.0*t44*t24-4.0*t27*t8*
                        t296*t30-2.0*t253*t317*t30-2.0*t65*t2*t23*t6*t30+2.0*t3*t23*t14*t30-2.0*t12*t19
                        *t331*t6+2.0*t335*t331*t30-2.0*t201*t339+2.0*t201*t342+2.0*t201*t345+2.0*t86*
                        t222;
                t354 = 1/(t95+t139+t167+t193+t219+t260+t303+t350);
                t361 = t22*t4+t29*t8+t296-t23*t4-t30*t8-h[5];
                t365 = t253*t18-t3*t23-t335*t2+t201+t339-t342-t345+t202;
                t374 = t66-2.0*t243+t49+t90-2.0*t225+t100+t35+t39+t41+t43+t45+t46;
                err[0] = sqrt((t17*t47*t354-t361*t365*t354)*t17+(-t17*t365*t354+t361*t374*
                                                                 t354)*t361);
                return;
        }
}
