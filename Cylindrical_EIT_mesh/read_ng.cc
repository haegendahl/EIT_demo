
#include "mex.h"
#include "readvolmesh.hh"

using namespace std;
mxArray* createNodeStruct(uint32 size);

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  char *ibuf;
  int buflen, order, dims[2];
  int isok = 0;

  /* check proper input and output */
  if(nrhs!=2)
    mexErrMsgTxt("2 inputs required.");
  else if(nlhs != 4)
    mexErrMsgTxt("4 output arguments required.");
  else if(!mxIsChar(prhs[0]))
    mexErrMsgTxt("Input 1 must be a string.");

  // get the length of the input string 
  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  

  // copy the string data from prhs[0] into a C string input_ buf.    
  ibuf = mxArrayToString(prhs[0]);
    
  if(ibuf == NULL) 
    mexErrMsgTxt("Could not convert input to string.");

    
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
      mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input 2 must be a scalar.");
  }
  
  order = (int) mxGetScalar(prhs[1]);  
  if (!((order == 1) | (order == 2))) {
    mexErrMsgTxt("Order must be 1 or 2.");
  }
 
  volmesh::ReadVolmesh netgenhila(order);

  netgenhila.open(ibuf);
  if (!netgenhila.isvolmesh()) {
    mexErrMsgTxt("Could not open file!");
  }
  
  isok = netgenhila.getMeshdata();
  if (!isok) mexErrMsgTxt("Could not read meshdata");

  // Volume element (H) output
  dims[1] = netgenhila.volsize(); dims[0] = 4 + 6*(order-1);
  plhs[0] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);

  // Node (g) output
  plhs[1] = mxCreateDoubleMatrix(3, netgenhila.nodesize(), mxREAL);
  
  // surface data (surf-ind, bcnr, h) output
  dims[1] = netgenhila.surfsize(); dims[0] = 3*order+2;
  plhs[2] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);

  // populate output matrices
  netgenhila.copyVolElements((uint32 *) mxGetData(plhs[0]));
  netgenhila.copySurfData((uint32 *) mxGetData(plhs[2]));
  netgenhila.copyNodes(mxGetPr(plhs[1]));

  plhs[3] = createNodeStruct(netgenhila.nodesize());
}

mxArray* createNodeStruct(uint32 size) {
  mxArray *node_struct;
  const char *field_names[] = {"Topology", "ElementConnection"};
  int dims[2];
  int i, t_field, e_field;
    
  dims[0] = (int) size;
  dims[1] = 1;

    /* Create a 1-by-n array of structs. */ 
  node_struct = mxCreateStructArray(2, dims, 2, field_names);

    /* This is redundant, but here for illustration.  Since we just
       created the structure and the field number indices are zero
       based, name_field will always be 0 and phone_field will always
       be 1 */
  t_field = mxGetFieldNumber(node_struct,"Topology");
  e_field = mxGetFieldNumber(node_struct,"ElementConnection");

    /* Populate the name and phone fields of the phonebook structure. */ 
  // mxArray *field_value;

  //field_value = mxCre   *mxGetPr(field_value) = friends[i].phone;

  for (i=0; i<4; i++) {
    mxSetFieldByNumber(node_struct,i,t_field,mxCreateString("hello"));
      mxSetFieldByNumber(node_struct,i,e_field,mxCreateDoubleMatrix(1,1,mxREAL));
  }
  
  return node_struct;
}
