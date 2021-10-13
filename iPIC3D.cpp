/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
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


#include "MPIdata.h"
#include "iPic3D.h"
#include "debug.h"
#include "TimeTasks.h"
#include <stdio.h>
#include <sys/time.h>

using namespace iPic3D;

int main(int argc, char **argv) {

 MPIdata::init(&argc, &argv);
 {
  iPic3D::c_Solver KCode;
  KCode.Init(argc, argv);

  timeTasks.resetCycle();
  KCode.CalculateMoments();

  // start measuring time
  //struct timeval begin, end, t0_in, t1_in;
  //gettimeofday(&begin, 0);

  for (int i = KCode.FirstCycle(); i < KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0)
      printf(" ======= Cycle %d ======= \n",i);
   
    //gettimeofday(&t0_in, 0);

    timeTasks.resetCycle();
    KCode.CalculateField(i);
    KCode.ParticlesMover(i);
    if (KCode.get_myrank() == 0)
      printf("** Particle Mover done **\n");
    KCode.CalculateB();
    KCode.CalculateMoments();
    KCode.WriteOutput(i);
    // print out total time for all tasks
    timeTasks.print_cycle_times(i);

    //gettimeofday(&t1_in, 0);

    /*if (KCode.get_myrank() == 0){
      long seconds = t1_in.tv_sec - t0_in.tv_sec;
      long microseconds = t1_in.tv_usec - t0_in.tv_usec;
      double elapsed = seconds + microseconds*1e-6;
      printf("Time rank%i measured: %.5f seconds.\n", KCode.get_myrank(), elapsed);
    }*/

  }

  // Stop measuring time and calculate the elapsed time
  //gettimeofday(&end, 0);
  //long seconds = end.tv_sec - begin.tv_sec;
  //long microseconds = end.tv_usec - begin.tv_usec;
  //double elapsed = seconds + microseconds*1e-6;
  //printf("Time rank%i measured: %.5f seconds.\n", KCode.get_myrank(), elapsed);

  KCode.Finalize();
 }
 // close MPI
 MPIdata::instance().finalize_mpi();

 return 0;
}
