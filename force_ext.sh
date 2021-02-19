gcc -shared -o Frc_ext.so -DUSE_TCL_STUBS -I$TCLINC -L$TCLLIB -ltclstub8.6 -fPIC frc_ext.c
