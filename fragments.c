#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main ()
{

  //************************************
  // input
  //************************************

  // define input/output file names
  char inputFile[]="fragment_id.txt";
  char outputFile[]="cluster_size.txt";
  char plotFile[]="plot.txt";
  char histFile[]="hist.txt";
  char numbFile[]="clus_nums.txt";

  // misc integers
  int i;

  // lines to ignore in read
  int  nlines=7;
  char **misclines;

  // output array size (should be geq # max # of fragments)
  int maxfrags=20000;
  int output[maxfrags][2];
  memset(output, 0, maxfrags*2*sizeof(int));

  misclines = (char **)malloc(nlines*sizeof(char *));
  for(i=0; i<nlines; i++)
    misclines[i] = (char *)malloc(256*sizeof(char));

  strcpy(misclines[0],"ITEM: BOX BOUNDS pp pp pp\n");
  strcpy(misclines[1],"0.0000000000000000e+00 6.0000000000000000e+01\n");
  strcpy(misclines[2],"0.0000000000000000e+00 4.2000000000000000e+01\n");
  strcpy(misclines[3],"0.0000000000000000e+00 2.0000000000000000e+01\n");
  strcpy(misclines[4],"0.0000000000000000e+00 5.0000000000000000e+01\n");
  strcpy(misclines[5],"ITEM: ATOMS id c_frag\n");
  strcpy(misclines[6],"ITEM: ATOMS id c_cluster\n");

  //************************************ 
  // working code
  // DO NOT CHANGE
  //************************************

  // open input file read
  FILE* fin=fopen(inputFile,"r");

  // define output file handles
  FILE* fout=fopen(outputFile,"a");
  FILE* fplot=fopen(plotFile,"a");
  FILE* fhist=fopen(histFile,"a");
  FILE* fnumb=fopen(numbFile,"a");

  int atoms_assigned=0;
  int atoms=0;
  int timestep=0;
  int count=2;    // initialized to two for matching at initial time step

  // read file line by line
  char* line=NULL;
  size_t len=120;
    

  while(getline(&line,&len,fin)!=-1)
  {
    // if TIMESTEP
    //   check total number of atoms, print output, reset output structures
    if(strcmp(line,"ITEM: TIMESTEP\n")==0)
    {
      // count should be atom plus two lines for line defining number of atoms
      if(count!=atoms+2)
      {
        printf("Atom count does not match\n");
        printf("Atoms=%d Count=%d\n",atoms,count);
        break;
      }
      // check if initial step
      if (count!=2)
        print_one(timestep,output,fout,fplot,fhist,fnumb);
      // reset output file
      memset(output, 0, sizeof(output[0][0])*maxfrags*2);
      // get new timestep
      getline(&line,&len,fin);
      timestep=atoi(line);
      printf("Timestep: %d\n",timestep);
      // reset count to zero
      count=0;
      continue;
    }

    // if misc lines continue
    if(ismiscline(line,misclines,nlines))
      continue;

    // if ATOM number
    //   if atoms unassigned, assign number of atoms
    //   continue
    if(strcmp(line,"ITEM: NUMBER OF ATOMS\n")==0)
    {
      getline(&line,&len,fin);
      count+=2;
      if(!atoms_assigned)
      {
        atoms=atoi(line);
        atoms_assigned=1;
      }
      continue;
    }

    // extract fragment id
    int tmp=0;
    int index=0;
    int frag_id;
    char* token = strtok(line," ");
    while(token!=NULL) 
    {
      if(tmp==0)
        index=atoi(token);
      if(tmp==1) 
        frag_id=atoi(token);

      tmp++;
      token=strtok(NULL," ");
    }
    if(tmp==2)
      count++;
    else
    {
      printf("Unexpected line: %s\n",line);
    }

    // add this fragment index to output file
    add_output(frag_id, output);

  }
  // Extra print for last step
  print_one(timestep,output,fout,fplot,fhist,fnumb);

  // free & fclose
  for(i=0; i<nlines; i++)
    free(misclines[i]);
  free(misclines);
  fclose(fin);
  fclose(fout);
  fclose(fplot);
  fclose(fhist);
  fclose(fnumb);

  return 1;
}

// is line one of the misc lines
int ismiscline(char* line, char** misclines, int nlines)
{
  int i;
  for(i=0; i<nlines; i++)
    if(strcmp(line,misclines[i])==0)
      return 1;
  
  return 0;
}

// print data from one step to outputfile
print_one(int timestep, int output[][2], FILE* fout, FILE* fplot, FILE* fhist, FILE* fnumb)
{
  // print file headers
  fprintf(fout,"Timestep: %d\n",timestep);
  fprintf(fhist,"Timetep: %d\n",timestep);
  fprintf(fnumb,"%d ",timestep);

  int i=0,j=0,k=0;
  int max=0;
  int max_one=0;
  int oldFragment=0;
  
  // max number of distinct fragment sizes
  int maxFrags=10000;

  // data strcuture for storing histogram at current time step
  unsigned long int** hist = (unsigned long int**)malloc(maxFrags * sizeof(unsigned long int*));
  for(j=0; j<maxFrags; j++)
  {
    hist[j]=(unsigned long int*)malloc(2 * sizeof(unsigned long int));
    memset(hist[j], 0, sizeof(unsigned long int) * 2);
  }

  // analyze fragment data
  while(output[i][0]!=0)
  {
    fprintf(fout,"%d %d\n",output[i][0],output[i][1]);
    if(output[i][1]>max)
    {
      max_one=max;
      max=output[i][1];
    }
    else if(output[i][1]>max_one)
    {
      max_one=output[i][1];
    }

    // add fragment size to histogram
    while(hist[k][0]!=0)
    {
      if(hist[k][0]==output[i][1])
      {
        oldFragment=1;
        break;
      }
      k++;
    }
    if(!oldFragment)
      hist[k][0]=output[i][1];
    else
      oldFragment=0;
    hist[k][1]++;
    k=0;


    i++;
  }

  unsigned long int DPn_num=0, DPw_num=0, DPz_num=0, DPn_r_num=0, DPw_r_num=0, DPz_r_num=0;
  unsigned long int totalFrags=0, n_total=0;

  // print analyzed data to files
  while(hist[k][0]!=0)
  {
    DPn_num+=hist[k][0]*hist[k][1];
    DPw_num+=hist[k][0]*hist[k][0]*hist[k][1];
    DPz_num+=hist[k][0]*hist[k][0]*hist[k][0]*hist[k][1];
    if(hist[k][0]!=max)
    {
      DPz_r_num+=hist[k][0]*hist[k][0]*hist[k][0]*hist[k][1];
      DPw_r_num+=hist[k][0]*hist[k][0]*hist[k][1];
      DPn_r_num+=hist[k][0]*hist[k][1];
    }
    // ignore only largest
    //   later there can be multiple frgments with size==largest size
    //   hence use hist -1 to ignore only largest
    else
    {
      DPz_r_num+=hist[k][0]*hist[k][0]*hist[k][0]*(hist[k][1]-1);
      DPw_r_num+=hist[k][0]*hist[k][0]*(hist[k][1]-1);
      DPn_r_num+=hist[k][0]*(hist[k][1]-1);
    }
    totalFrags+=hist[k][1];

    fprintf(fhist,"%d %d\n", hist[k][0], hist[k][1]);
    n_total+=hist[k][0]*hist[k][1];
    k++;
  }
  fprintf(fhist,"Total atoms: %d\n", n_total);
  fprintf(fhist,"Total mols : %d\n", totalFrags);
  fprintf(fnumb,"%d\n", totalFrags);

  double dpn  =((double)DPn_num)/((double)totalFrags);
  double dpw  =((double)DPw_num)/((double)DPn_num);
  double dpz  =((double)DPz_num)/((double)DPw_num);
  double dpw_r=((double)DPw_r_num)/((double)DPn_r_num);
  double dpz_r=((double)DPz_r_num)/((double)DPw_r_num);
  double pdi  =dpw/dpn;
  double q    =0.5*(dpw-1);

  fprintf(fplot,"%lu %f %f %f %f %d %f %f %f\n",timestep, dpn, dpw, dpw_r, pdi, max, q, dpz, dpz_r);

  // close output files and dynamic alloc
  for(j=0; j<maxFrags; j++)
    free(hist[j]);
  free(hist);
}

// add frag id to output array
add_output(int frag_id, int output[][2])
{
  int i=0;
  int isold=0;
  while(output[i][0]!=0)
  {
    if(frag_id==output[i][0])
    {
      isold=1;
      break;
    }
    i++;
  }

  if(isold)
    output[i][1]++;
  else
  {
    output[i][0]=frag_id;
    output[i][1]++;
  }
}
