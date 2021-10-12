#cd source_code
#g++ efficiency.C Controll.C makeHists.C efficiency.C -g -o efficiency.exe `root-config --cflags --glibs`
g++ -g -std=c++17 Controll.C makeHists.C efficiency.C -I$ROOTSYS/include `root-config --libs ` -lMinuit -lGenVector -o efficiency.exe
#cd ../run/4l
