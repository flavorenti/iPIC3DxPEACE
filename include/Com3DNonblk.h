/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * iPIC3DxPEACE was later developed by Federico Lavorenti for planet applications.
 * Publications that use results from iPIC3DxPEACE need to properly cite these works: 
 * (1) 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 * (2) ' F. Lavorenti, P. Henri ... '
 *
 * Copyright 2022 Observatoire de la Cote d'Azur, Nice, France
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef Com3DNonblk_H
#define Com3DNonblk_H

#include "arraysfwd.h"
#include "ipicfwd.h"
#include "mpi.h"
#include "BcFields3D.h"
#include "VCtopology3D.h"
#include "ipicdefs.h"
#include "Alloc.h"
#include "debug.h"
#include "EMfields3D.h"

/*! @file 
 *  Library for Nonblocking Halo Exchange.
 *  Developers: Stefano Markidis, Ivy Bo Peng (Feb 2015)
 */

/*! @brief Non-Blocking (NB) derived Halo values communication among neighbour MPI procs.
 *
 * ... missing a long description ...
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3-pointer providing values to be communicated 
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 * @param[in] isCenterFlag True for mesh centers and False for mesh Nodes 
 * @param[in] isFaceOnlyFlag True when using BoxStencil and False otherwise 
 * @param[in] needInterp True to use interpolation and False otherwise
 * @param[in] isParticle True when using pcls communicator and False when using fields comm.
 */
void NBDerivedHaloComm(int nx, int ny, int nz, 
		       double ***vector,
		       const VirtualTopology3D * vct, 
		       EMfields3D *EMf,
		       bool isCenterFlag, 
		       bool isFaceOnlyFlag, 
		       bool needInterp, 
		       bool isParticle);

/*! @brief Communicate values in ghost cells among MPI procs (centers of the mesh).
 *
 * Function make communications between neighbours MPI procs in the 3D cartesian topology. 
 * The values at the centers of the grid mesh in the ghost cells are exchanged among neighbour MPI procs.
 * Centers meaning values like V[i+1/2]. Nodes meaning values V[i].
 * Inside this function there are calls to two other functions: NBderivedHaloComm (that make the MPI comm) and BCface (that apply BC).
 * The behaviour of NBDerivedHaloComm is specified by 4 flags passed as input parameters. In this case TFFF.
 * BCface is used to ensure good-beahviour of box-boundary MPI procs.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateCenterBC(int nx, int ny, int nz, 
		         arr3_double vector, 
			 int bcFaceXright, int bcFaceXleft, 
			 int bcFaceYright, int bcFaceYleft, 
			 int bcFaceZright, int bcFaceZleft, 
			 const VirtualTopology3D * vct, 
			 EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (centers of the mesh) for pcls.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case TFFT.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateCenterBC_P(int nx, int ny, int nz, 
		           arr3_double vector, 
			   int bcFaceXright, int bcFaceXleft, 
			   int bcFaceYright, int bcFaceYleft, 
			   int bcFaceZright, int bcFaceZleft, 
			   const VirtualTopology3D * vct, 
			   EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (centers of the mesh) using stencil.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case TTFF.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateCenterBoxStencilBC(int nx, int ny, int nz, 
		                   arr3_double vector, 
				   int bcFaceXright, int bcFaceXleft, 
				   int bcFaceYright, int bcFaceYleft, 
				   int bcFaceZright, int bcFaceZleft, 
				   const VirtualTopology3D * vct, 
				   EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (nodes of the mesh) using stencil for pcls.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case FTFT.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, 
		                   arr3_double vector, 
				   int bcFaceXright, int bcFaceXleft, 
				   int bcFaceYright, int bcFaceYleft, 
				   int bcFaceZright, int bcFaceZleft, 
				   const VirtualTopology3D * vct, 
				   EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs with interpolation.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case TFTT. No BC in this case.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3-pointer that provides the values for communication 
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateInterp(int nx, int ny, int nz, 
		       double*** vector, 
		       const VirtualTopology3D * vct, 
		       EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (nodes of the mesh) for pcls.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case FFFT. No BC in this case.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3-pointer that provides the values for communication 
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateNode_P(int nx, int ny, int nz, 
		       double*** vector, 
		       const VirtualTopology3D * vct, 
		       EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (nodes of the mesh).
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case FFFF.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateNodeBC(int nx, int ny, int nz,
	               arr3_double vector, 
		       int bcFaceXright, int bcFaceXleft, 
		       int bcFaceYright, int bcFaceYleft, 
		       int bcFaceZright, int bcFaceZleft, 
		       const VirtualTopology3D * vct, 
		       EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (nodes of the mesh) for pcls.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case FFFT.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateNodeBC_P(int nx, int ny, int nz, 
		         arr3_double vector, 
			 int bcFaceXright, int bcFaceXleft, 
			 int bcFaceYright, int bcFaceYleft, 
			 int bcFaceZright, int bcFaceZleft, 
			 const VirtualTopology3D * vct, 
			 EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (nodes of the mesh) using stencil.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case FTFF.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateNodeBoxStencilBC(int nx, int ny, int nz, 
		                 arr3_double vector, 
				 int bcFaceXright, int bcFaceXleft, 
				 int bcFaceYright, int bcFaceYleft, 
				 int bcFaceZright, int bcFaceZleft, 
				 const VirtualTopology3D * vct, 
				 EMfields3D *EMf);

/*! @brief Communicate values in ghost cells among MPI procs (centers of the mesh) using stencil for pcls.
 *
 * Function similar to communicateCenterBC. 
 * Only difference is the flag passed to NBDerivedComm. In this case TTFT.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] vector 3D built-in type, provide values for communcation 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] EMf an instance of the electromagnetic field class
 */
void communicateCenterBoxStencilBC_P(int nx, int ny, int nz, 
		                     arr3_double vector, 
				     int bcFaceXright, int bcFaceXleft, 
				     int bcFaceYright, int bcFaceYleft, 
				     int bcFaceZright, int bcFaceZleft, 
				     const VirtualTopology3D * vct, 
				     EMfields3D *EMf);
#endif
