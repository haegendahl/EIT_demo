#ifndef READVOLMESH
#define READVOLMESH

// default size of read buffer

//#include "bufreader.hh"
#include <fstream>
#include <iostream>
#include "volmeshtypes.hh"

namespace volmesh {

class ReadVolmesh {
public:
  // ReadVolmesh();  // lienee turhake
  ReadVolmesh(char baseorder);
  ~ReadVolmesh();
  
  int open(char *fname);
  
  char getFilename();
  
  int getMeshdata();
  int getHeader();
  int getSurfaceElements();
  int getVolumeElements();
  int getNodes();
  int getDomains();
  int getBoundaries();
  
  uint32 volind();
  uint32 volind(uint32 elid, byte np);
  uint32 volsize();
  uint32 volelements();

  uint32 surfsize();
  uint32 nodesize();

  bool isvolmesh();

  void close();
  void print_surfind();
  void print_bcnr();
  void print_surf(); 
  void print_volume(); 
  void print_nodes(); 

  void copyVolElements(uint32 *mem);
  void copySurfData(uint32 *mem);
  void copyNodes(double *mem);

  void elementConnection();

  template <class T>
  T getElement(T eind) {
    if (eind<1) eind=1;
    if (eind>nH) eind=nH;
    return H[eind];
  }

private:
  std::ifstream file;
  char *filename, *buffer;
  uint32 *H, *h, *surf, *bcnr;
  double *g; 

  uint32 nH, nh, ng, epoints, volpos, sizeH, sizeh, sizeg;
  char order;

  bool ismesh, issurface, issurfind, isvolume, isnode, isfail;

  int findHeader(char* str);
  char getchar(char *str, int slen);

  void readelements(uint32 *data, uint32 Ndata, uint32 np, int hpos, int offset);
  void readelements(double *data, uint32 Ndata, uint32 np);

};
 


}
#endif
