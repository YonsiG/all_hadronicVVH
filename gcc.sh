g++ -g -std=c++17 Controll.C makeHists.C efficiency.C -I$ROOTSYS/include `root-config --libs ` -lMinuit -lGenVector -o efficiency.exe
