
/*  step_CNT
 * VILLIN
 * Tcl 8.x wrapper for the Step function.
 * 1/1+exp(beta(dist-landa(d0)))
 */
 #include <math.h>
 #include <tcl.h>
 #include <stdlib.h>
 #include <stdio.h>
 

 static double  Q_fun, state[156], vecsub[3]={0}, grad[3]={0}; 
 static double  Kw = 1400.0, wr, vecscl, Fcons, d_d0, dist, denom, Exp, beta = 5.0, sysiz = 156.0, sysiz_invrt = 0.006410256410256;
 static int     pair[156][2], sortpair[156][2], usedatoms[92]; /* 156 is the number of pairs, and 92 is the number of atoms participat in the simulation */
 static int     stop = 1, pull; 
 static int     simNum,ts,oldTS = -1; /* ts should be equal to the firsttimestep in the config file */ 
 static FILE    *gFile; 
 static char    buffer1 [20], buffer2 [20], buffer3 [100], buffer5 [10];
 
 static int     replica_id;

 frc_CNT(ClientData clientData, Tcl_Interp *interp,
        int objc, Tcl_Obj *CONST objv[]) {
    
    double  data1, data2;
    double  plst[92][3]={0};
    double  N = 0.0; 
    int     error, num1, num2, i;
    

    Tcl_Obj *value1 = Tcl_GetVar2Ex(interp, "wr", NULL, TCL_GLOBAL_ONLY);
    Tcl_GetDoubleFromObj(interp,value1,&wr);
    Tcl_Obj *value3 = Tcl_GetVar2Ex(interp, "ts", NULL, TCL_GLOBAL_ONLY);
    Tcl_GetIntFromObj(interp,value3,&ts);
  
    if (stop == 1) {
      Tcl_Obj *value2 = Tcl_GetVar2Ex(interp, "replica_id", NULL, TCL_GLOBAL_ONLY);
      Tcl_GetIntFromObj(interp,value2,&replica_id);
      Tcl_Obj *value4 = Tcl_GetVar2Ex(interp, "i_job", NULL, TCL_GLOBAL_ONLY);
      Tcl_GetIntFromObj(interp,value4,&simNum);
      
        Tables();
        stop = 0;
    }
        
    Tcl_Obj **list2, **list1;
    Tcl_Obj *Coord1, *Coord2;
     
    for (i=0 ; i < 156 ; i++) { 
        num1 = 0; num2 = 0;
        snprintf ( buffer1, 20, "atmcrd(%d)" , pair[i][0] );
        snprintf ( buffer2, 20, "atmcrd(%d)" , pair[i][1] );
        
        Coord1 = Tcl_GetVar2Ex(interp, buffer1 , NULL, TCL_GLOBAL_ONLY);
        Coord2 = Tcl_GetVar2Ex(interp, buffer2 , NULL, TCL_GLOBAL_ONLY);
        if (Tcl_ListObjGetElements(interp, Coord1, &num1, &list1) != TCL_OK) return TCL_ERROR;
        if (Tcl_ListObjGetElements(interp, Coord2, &num2, &list2) != TCL_OK) return TCL_ERROR;
        
        if (Tcl_GetDoubleFromObj(interp, list1[0], &data1) != TCL_OK) return TCL_ERROR;
        if (Tcl_GetDoubleFromObj(interp, list2[0], &data2) != TCL_OK) return TCL_ERROR;
        vecsub[0] = data1 - data2 ;
        if (Tcl_GetDoubleFromObj(interp, list1[1], &data1) != TCL_OK) return TCL_ERROR;
        if (Tcl_GetDoubleFromObj(interp, list2[1], &data2) != TCL_OK) return TCL_ERROR;
        vecsub[1] = data1 - data2 ;
        if (Tcl_GetDoubleFromObj(interp, list1[2], &data1) != TCL_OK) return TCL_ERROR;
        if (Tcl_GetDoubleFromObj(interp, list2[2], &data2) != TCL_OK) return TCL_ERROR;
        vecsub[2] = data1 - data2 ;
        
        dist = sqrt(vecsub[0]*vecsub[0]+vecsub[1]*vecsub[1]+vecsub[2]*vecsub[2]);
        d_d0 = dist - state[i];
	
	if  (d_d0 <= 3.5 && d_d0 >= -3.5) {
	   Exp    = exp(beta * d_d0);
           denom  = 1 + Exp;
           Q_fun  = 1 / (sysiz * denom);
           N      = Q_fun + N;
	   vecscl = -beta * Exp * Q_fun / denom / dist;
	}
	else if (d_d0 < -3.5) {
	   N  = sysiz_invrt + N;
	   continue;
	}
	else  {
	   continue;
	}

        grad[0]  = vecsub[0] * vecscl;
        grad[1]  = vecsub[1] * vecscl;
        grad[2]  = vecsub[2] * vecscl;
        
        plst[sortpair[i][0]][0] = plst[sortpair[i][0]][0] + grad[0];
        plst[sortpair[i][0]][1] = plst[sortpair[i][0]][1] + grad[1];
        plst[sortpair[i][0]][2] = plst[sortpair[i][0]][2] + grad[2];        
        plst[sortpair[i][1]][0] = plst[sortpair[i][1]][0] - grad[0];
        plst[sortpair[i][1]][1] = plst[sortpair[i][1]][1] - grad[1];
        plst[sortpair[i][1]][2] = plst[sortpair[i][1]][2] - grad[2];
    }
    
    
    if (oldTS < ts){
      fprintf (gFile,"%0.5f\n",N);
      oldTS = ts;
    }
    snprintf ( buffer5, 10, "%f" ,N);
    Tcl_SetVar(interp, "n", buffer5, TCL_GLOBAL_ONLY);

    
    if (ts % 100000 == 0) {
        fflush (gFile);
	if (fflush(gFile)) { 
	  perror ("Flush Error");
	  CkExit();
	}
    }
    Fcons = (wr - N) * Kw;
    ts++;
      
    for (i=0; i<92; i++) {
        snprintf ( buffer3, 100, "addforce %d {%.15f %.15f %.15f} " ,
                usedatoms[i], plst[i][0] * Fcons, plst[i][1] * Fcons, plst[i][2] * Fcons);
        Tcl_GlobalEval(interp, buffer3);
    }
        
 }





 int Frc_ext_Init(Tcl_Interp *interp) {
     Tcl_CreateObjCommand(interp, "frc", frc_CNT,
             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
     return TCL_OK;
 }








                     // \\
                    // * \\ 
                   // *** \\ 
                  // ***** \\ 
                 // ******* \\ 
                // ********* \\ 
               // *********** \\ 
              // ************* \\ 
             // *** Tables **** \\ 
            // ***************** \\ 
           // ******************* \\ 
          // ********************* \\ 
         // *********************** \\   
  // *************************************** \\     
 // **** Coordinate of Cristal Structure **** \\   
// ******************************************* \\ 

 static char    buffer4 [30];    
 int Tables() {
    
    FILE *pFile;
    pFile = fopen("input/state.dat","r");
    
    int i = 0; double num1;
    while (!feof(pFile)) {
        fscanf(pFile,"%lf",&num1);
        state[i] = num1;
        i++;
    }
    fclose(pFile);
    
// ******************************************* \\
// *************** Pair Atoms **************** \\
// ******************************************* \\
    
    pFile = fopen("input/pair.dat","r");
    
    i = 0; int num2; int k = 0;
    while (!feof(pFile)) {
        fscanf(pFile,"%d",&num2);
        if (k == 0)  {
            k = 1;
            pair[i][0] = num2;
        }
        else  {
            pair[i][1] = num2;
            k = 0;
            i++;
        }
    }
    fclose(pFile);
    
// ******************************************* \\
// ***** Ascending of pair atoms' Index ****** \\
// ******************************************* \\
    
    pFile = fopen("input/sort_pair.dat","r");
    
    i = 0; int num3; k = 0;
    while (!feof(pFile)) {
        fscanf(pFile,"%d",&num3);
        if (k == 0)  {
            k = 1;
            sortpair[i][0] = num3;
        }
        else  {
            sortpair[i][1] = num3;
            k = 0;
            i++;
        }
    }
    fclose(pFile);
    
// ******************************************* \\
// ******** used atoms in simulation ********* \\
// ******************************************* \\
    
    pFile = fopen("input/usedatoms.dat","r");
    
    i = 0; int num6;
    while (!feof(pFile)) {
        fscanf(pFile,"%d",&num6);
        usedatoms[i] = num6;
        i++;
    }
    fclose(pFile);
    
    
// ******************************************* \\
// ************** SumN.dat File ************** \\
// ******************************************* \\
    
    
    snprintf ( buffer4, 30, "output_%03d/sumN_%03d.%03d.dat" ,simNum ,simNum, replica_id );
    gFile = fopen (buffer4,"w");
 }
                        /*
  // *************************************** \\     
 // ****************** End ****************** \\   
// ******************************************* \\
          \\ *********************** //
           \\ ********************* // 
            \\ ******************* // 
             \\ ***************** // 
              \\ *************** //
               \\ ************* //
                \\ *********** //
                 \\ ********* //
                  \\ ******* // 
                   \\ ***** //
                    \\ *** //
                        */
			
			
			


