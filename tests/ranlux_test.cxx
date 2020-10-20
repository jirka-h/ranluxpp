/*************************************************************************
 * Copyright (C) 2018,  Alexei Sibidanov                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * This program shows usage of the RANLUX++ random number generator and  *
 * also benchmarks it.                                                   *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU Lesser General Public License as        *
 * published by the Free Software Foundation, either version 3 of the    *
 * License, or (at your option) any later version.                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 * Lesser General Public License for more details.                       *
 *                                                                       *
 * You should have received a copy of the GNU Lesser General Public      *
 * License along with this program.  If not, see                         *
 * <http://www.gnu.org/licenses/>.                                       *
 *************************************************************************/

#include "ranlux.h"
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <inttypes.h>

// test the performance by summing up 10^9 float numbers
// note the sum is always 2^24=16777216 due to rounding errors
template<class T>
void speedtest(){
  T g1(3124);
  int N = 1000*1000*1000;
  float sum = 0;
  printf("Summing up %d random floats...\n",N);
  for(int i=0;i<N;i++) sum += g1();
  printf("sum=%lf\n", (double)(sum));
}

// test the performance by skipping 10^9 states without actual
// delivery of random numbers to the user
template<class T>
void speedtest_nextstate(){
  T g1(3124);
  int N = 1000*1000, M = 1000;
  printf("Skipping %d states...\n",N*M);
  for(int i=0;i<N;i++) g1.nextstate(M);
  printf("Done.\n");
}

void test_sameseed(){
  int N = 1000*1000, M = 100;
  printf("Seed all SIMD generators with the same seed -- has to duplicate the scalar version.\n\n");

  ranluxI_scalar g1(3124);
  printf("Skipping %d states (scalar) ...\n",N*M);
  for(int i=0;i<N;i++) g1.nextstate(M);
  printf("Done.\n");
  for(int i=0;i<24;i++) printf("%f ", g1()); printf("\n\n");

  ranluxI_SSE g2(3124); g2.init(3124,1);
  printf("Skipping %d states (SSE2) ...\n",N*M);
  for(int i=0;i<N;i++) g2.nextstate(M);
  printf("Done.\n");
  for(int i=0;i<4*24;i++) printf("%f ", g2()); printf("\n\n");

#ifdef __AVX2__
  ranluxI_AVX g3(3124); g3.init(3124,1);
  printf("Skipping %d states (AVX2) ...\n",N*M);
  for(int i=0;i<N;i++) g3.nextstate(M);
  printf("Done.\n");
  for(int i=0;i<8*24;i++) printf("%f ", g3()); printf("\n\n");
#endif
}

// compare results with the original FORTRAN code:
// http://www.cpc.cs.qub.ac.uk/summaries/ACPR_v1_0.html or
// http://www.cpc.cs.qub.ac.uk/summaries/ACPR_v2_0.html
template<class ranlux_imp>
void original_test(){
  float rvec[1000];
  ranlux_imp a;

  auto print = [&](const char *title, const char *luxury){
    printf("  %s\n",title);
    a.ranlux(rvec,100); printf(" RANLUX %s   1-  5:\n",luxury);
    for(int i=0;i<5;i++) printf("%10.8f ",rvec[i]); printf("\n");
    a.ranlux(rvec,100); printf(" RANLUX %s 101-105:\n",luxury);
    for(int i=0;i<5;i++) printf("%10.8f ",rvec[i]); printf("\n");
  };

  print(" CALL RANLUX(RVEC,100)","default numbers");
  a.rluxgo(0,0,0,0);
  print(" CALL RLUXGO(0,0,0,0)","luxury level 0,");
  a.rluxgo(389,1,0,0);
  print(" CALL RLUXGO(389,1,0,0)","luxury p=389,");
  a.rluxgo(75,0,0,0);
  print(" CALL RLUXGO(75,0,0,0)","luxury p= 75,");

  int isdext[25];
  printf("  test restarting from the full vector\n");
  a.rluxut(isdext);
  printf("current RANLUX status saved:\n");
  for(int i=0;i<25;i++) {
    printf("%12d", isdext[i]); if(!((i+1)%5)) printf("\n");
  }
  printf("\n");
  print("","numbers");
  printf("   previous RANLUX status will be restored\n");
  a.rluxin(isdext);
  print("","numbers");

  printf("     test the restarting by skipping\n");
  a.rluxgo(4,7674985,0,0);
  int i1,i2,i3,i4;
  a.rluxat(i1,i2,i3,i4);
  printf("  RLUXAT values = %d %d %d %d\n",i1,i2,i3,i4);
  for(int i=0;i<10;i++) a.ranlux(rvec,1000);
  a.rluxat(i1,i2,i3,i4);
  printf("  RLUXAT values = %d %d %d %d\n",i1,i2,i3,i4);
  a.ranlux(rvec,200);
  printf("  Next and 200th numbers are: %10.6f %10.6f\n",rvec[0],rvec[199]);
  a.rluxgo(i1,i2,i3,i4);
  printf("  Next and 200th numbers are: %10.6f %10.6f\n",rvec[0],rvec[199]);
}

void output_to_file(const char * filename) {

  signal(SIGPIPE, SIG_IGN); 

  FILE * stream;
  stream = fopen (filename,"w");
  if ( !stream ) {
    perror("Error on fopen");
    return;
  }

  ranluxI_scalar g1(3124);
  const int steps = 2048;
  const double giga=1073741824;
//from 24 24bits integers we can get 24*24/32=18 32-bits integers
  const size_t word_size = 18;
  uint32_t *state = new uint32_t[word_size];
  const size_t N = word_size * steps;
  uint32_t *buf = new uint32_t[N];

  uint32_t *p = buf;
  size_t rc;
  uint64_t total=0;

  for(;;) {
    p=buf;
    for (int i=0;i<steps;++i) {
      g1.nextstate_and_get_uint32_vector(state);
      memcpy(p, state, word_size*sizeof(uint32_t));
      p+=word_size;
    }
    rc = fwrite(buf, sizeof(uint32_t), N, stream);
    total += rc;
    if ( rc < N ) {
      perror("fwrite");
      fprintf(stderr, "ERROR: fwrite - bytes written %zu, bytes to write %zu\n",
            rc * sizeof(uint64_t), N * sizeof(uint64_t));
      fprintf(stderr, "Total bytes written %" PRIu64 ", %g GiB\n", total*sizeof(uint64_t),(double)(total*sizeof(uint64_t))/giga );
      delete[] buf;
      return;
    }
  }
}

void usage(int argc, char **argv){
  printf("Program to test the performance of the optimized RANLUX implementations (with skipping).\n");
  printf("Usage: %s ntest\n",argv[0]);
  printf("  ntest: 0 -- perform self consistency test\n");
  printf("              (random numbers are the same as in the original FORTRAN code)\n");
  printf("         1 -- sum of 10^9 float random numbers with the scalar skipping\n");
  printf("         2 -- sum of 10^9 float random numbers with the SSE2 skipping\n");
  printf("         3 -- sum of 10^9 float random numbers with the AVX2 skipping\n");
  printf("         4 -- skip 10^9 states or 1*24*10^9 numbers with the scalar skipping\n");
  printf("         5 -- skip 10^9 states or 4*24*10^9 numbers with the SSE2 skipping\n");
  printf("         6 -- skip 10^9 states or 8*24*10^9 numbers with the AVX2 skipping\n");
  printf("         7 -- same seed for SIMD generators (consistency check)\n");
  printf("         8 -- perform self consistency test using LCG as a skipping engine\n");
  printf("              (random numbers are the same as in the original FORTRAN code)\n");
  printf("         9 -- output stream of 64-bit random numbers. Filename required.\n");
  printf("              Example: %s 6 >(PractRand-RNG_test stdin32 -tlmax 32T -multithreaded)\n", argv[0]);
}

int main(int argc, char **argv){
  if(argc==1||argc>3) { usage(argc,argv); return 0;}

  int ntest = atoi(argv[1]);
  if ( ntest == 0 ){
    original_test<ranluxI_James>();
  } else if(ntest == 1){
    speedtest<ranluxI_scalar>();
  } else if(ntest == 2){
    speedtest<ranluxI_SSE>();
  } else if(ntest == 3){
#ifdef __AVX2__
    speedtest<ranluxI_AVX>();
#endif
  } else if(ntest == 4){
    speedtest_nextstate<ranluxI_scalar>();
  } else if(ntest == 5){
    speedtest_nextstate<ranluxI_SSE>();
  } else if(ntest == 6){
#ifdef __AVX2__
    speedtest_nextstate<ranluxI_AVX>();
#endif
  } else if(ntest == 7){
    test_sameseed();
  } else if(ntest == 8){
    original_test<ranluxpp_James>();
  } else if(ntest == 9){
    if (argc !=3)  { usage(argc,argv); return 0;}
    output_to_file(argv[2]);
  } else {
    usage(argc,argv);
  }
  return 0;
}
