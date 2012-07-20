/*
 *  main.c
 *  
 *
 *  Created by Earl Lawrence on 9/17/09.
 *
 *  This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL) 
 *  under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program 
 *  are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to 
 *  copy and use this software without charge, provided that this Notice and any statement of authorship are 
 *  reproduced on all copies.  Neither the U.S. Government nor LANS makes any warranty, express or implied, 
 *  or assumes any liability or responsibility for the use of this software.
 *
 */
 
 
#include "stdio.h"
#include "math.h"
#include "stdlib.h"

main(int argc, char **argv) {
    int i,j,type=1, writeout=0;
    double xstar[6], ystar[2*1995], stuff[4];
    FILE *fp;

    char fname[256];

    // xstar contains the 5 emulator parameters plus the red shift.

    // M000 
    /*h_CMB
    if ((fp = fopen("xstar.dat","r"))==NULL) {
        printf("Cannot find inputs.\n");
        exit(1);
    }
    for(j=0; j<5; j++) {
        fscanf(fp, "%lf", &xstar[j]);
    }
    fclose(fp);
    */


    xstar[0] = 1.344000e-01;
    xstar[1] = 2.246000e-02;
    xstar[2] = 9.610000e-01;
    xstar[3] = 8.070000e-01;
    xstar[4] = -1.000;
    xstar[5] = 0.0;

  sprintf(fname,"pdelta_11");

//     printf("Enter filename for output: ");
//     scanf("%s",fname);
//     printf("output will be written to: %s\n",fname);
// 
//     printf("Enter (Omega_m)*h^2: ");
//     scanf("%lf",&xstar[0]);
//     printf("%g\n",xstar[0]);
//     printf("Enter (Omega_b)*h^2: ");
//     scanf("%lf",&xstar[1]);
//     printf("%g\n",xstar[1]);
//     printf("Enter n_s: ");
//     scanf("%lf",&xstar[2]);
//     printf("%g\n",xstar[2]);
//     printf("Enter sigma_8: ");
//     scanf("%lf",&xstar[3]);
//     printf("%g\n",xstar[3]);
//     printf("Enter w: ");
//     scanf("%lf",&xstar[4]);
//     printf("%g\n",xstar[4]);
//     printf("Enter z: ");
//     scanf("%lf",&xstar[5]);
//     printf("%g\n",xstar[5]);
//     
//     printf("Enter output type (0: Delta^2/k^1.5; 1: Delta^2; 2: P(k)): ");
//     scanf("%li",&type);
//     printf("%i\n", type);
//     
type=1;
    getH0fromCMB(xstar, stuff);
    emu(xstar, ystar, &type);

    if ((fp = fopen(fname,"w"))==NULL) {
        printf("cannot open %s \n",fname);
        exit(1);
    }
    
    fprintf(fp, "# Parameters:\n");
    fprintf(fp, "# (Omega_m)*h^2 = %f, (Omega_b)*h^2 = %f, n_s = %f, sigma_8 = %f, w = %f\n", xstar[0], xstar[1], xstar[2], xstar[3], xstar[4]);
    fprintf(fp, "# z = %f\n", xstar[5]);
    fprintf(fp, "#\n");
    fprintf(fp, "# sound horizon: %g Mpc\n", stuff[0]);
    fprintf(fp, "# z_lss: %g\n", stuff[1]);
    fprintf(fp, "# d_lss from CMB: %g Mpc\n", stuff[2]);
    fprintf(fp, "# dimensionless Hubble parameter derived from CMB constraint: %1.3f\n",stuff[3]);
        
    switch(type) {
        default:
            fprintf(fp, "#\n# k[1/Mpc] Delta^2 / k^1.5:\n");
            break;
        case 1:
            fprintf(fp, "#\n# k[1/Mpc] Delta^2:\n");
            break;
        case 2:
            fprintf(fp, "#\n# k[1/Mpc] P(k):\n");
            break;
    }
    for(j=0; j<1995; j++) {
        for(i=0; i<=1; i++) {
            fprintf(fp, "%f ", ystar[i*1995+j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}
