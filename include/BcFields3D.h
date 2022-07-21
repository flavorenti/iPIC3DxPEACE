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

#ifndef BcFields_H
#define BcFields_H

#include "VCtopology3D.h"

/*! @file 
 *  Library to manage box boundary conditions.
 *  Developers: Stefano Markidis, Giovanni Lapenta (Jan 2009)
 */

/*! @brief Set the boundary condition on 3D box last boundary cell using field communicator.
 *
 * Function defining the rule to give value to last cell in your domain (e.g. V[0][j][k] in the case of Xleft).
 * This rule can be (i) asymmetric V[0]=-V[1] using flag=0, (ii) null V[0]=0 using flag=1, (iii) symmetric V[0]=V[1] 
 * using flag=2. The value of flag must be specified in the input file.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[out] vector 3-pointer to be updates at boundaries 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 */
void BCface(int nx, int ny, int nz, 
	    double ***vector, 
	    int bcFaceXright, int bcFaceXleft, 
	    int bcFaceYright, int bcFaceYleft, 
	    int bcFaceZright, int bcFaceZleft, 
	    const VirtualTopology3D * vct);

/*! @brief Set the boundary condition on 3D box last boundary cell using pcls communicator.
 *
 * Same as BCface but uses pcls communicator vct->getXright_neighbor_P() instead of the one for
 * fields vct->getXright_neighbor().
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[out] vector 3-pointer to be updates at boundaries 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology 
 */
void BCface_P(int nx, int ny, int nz, 
	      double ***vector, 
	      int bcFaceXright, int bcFaceXleft, 
	      int bcFaceYright, int bcFaceYleft, 
	      int bcFaceZright, int bcFaceZleft, 
	      const VirtualTopology3D * vct);

/*! @brief set the boundary condition on 4D box last boundary cell using field communicator.
 *
 * Same as BCface but acts on a 4D array (nx,ny,nz,ns), also one dimension
 * for number of species.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] ns number of the species from inp file
 * @param[out] vector 3-pointer to be updates at boundaries 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology 
 */
void BCface(int nx, int ny, int nz, int ns, 
	    double ****vector, 
	    int bcFaceXright, int bcFaceXleft, 
	    int bcFaceYright, int bcFaceYleft, 
	    int bcFaceZright, int bcFaceZleft, 
	    const VirtualTopology3D * vct);

/*! @brief set the boundary condition on 4D box last boundary cell using pcls communicator.
 *
 * Same as BCface_P but acts on a 4D array (nx,ny,nz,ns), also one dimension
 * for number of species.
 *
 * @param[in] nx dimension of the box in X from inp file
 * @param[in] ny dimension of the box in Y
 * @param[in] nz dimension of the box in Z
 * @param[in] ns number of the species from inp file
 * @param[out] vector 3-pointer to be updates at boundaries 
 * @param[in] bcFaceXright flag for BC anti-sunward (outflow side) from inp file
 * @param[in] bcFaceXleft flag for BC sunward (inflow side)
 * @param[in] bcFaceYright flag for BC (dusk side)
 * @param[in] bcFaceYleft flag for BC (dawn side)
 * @param[in] bcFaceZright flag for BC (north side)
 * @param[in] bcFaceZleft flag for BC (south side)
 * @param[in] vct built-in type defining the 3D MPI cartesian topology 
 */
void BCface_P(int nx, int ny, int nz, int ns, 
	      double ****vector, 
	      int bcFaceXright, int bcFaceXleft, 
	      int bcFaceYright, int bcFaceYleft, 
	      int bcFaceZright, int bcFaceZleft, 
	      const VirtualTopology3D * vct);

#endif

