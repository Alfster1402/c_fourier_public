#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979

//////////////////////// Defining the functions h_1 and h_2 /////////////////////////////////////

double real_h_1(int k) 
{
  return cos(k*PI/50) + cos(k*PI/10); // replacing t with value that corresponds to k, the value 
    // being summed over, which is dictated by the number of sample points N = 100
}

double imag_h_1(int k) 
{
  return sin(k*PI/50) + sin(k*PI/10);
}

double h_2(int k)
{
  return exp(((k*PI/50 - PI)*(k*PI/50 - PI))/2);
}

///////////// Defining the complex values of h_3 and the transformed values H(omega) ////////////

typedef struct h_3
{
  double real_h_3, imag_h_3;
} h_3;

typedef struct H_1
{
  double real_H_1, imag_H_1;
} H_1;

typedef struct H_2
{
  double real_H_2, imag_H_2;
} H_2;

typedef struct H_3 
{
  double real_H_3, imag_H_3;
} H_3;
  
////////////////// Defining the discrete Fourier transforms using structs ///////////////////////

struct H_1 DFT_h_1(struct H_1 result, int n) // Discrete Fourier transform for h_1
{ 
  int k = 0;
  result.real_H_1 = 0;
  result.imag_H_1 = 0;
  while (k < 100)
  {
    result.real_H_1 = result.real_H_1 + cos(k*PI*(1-n)/50) + cos(k*PI*(5-n)/50);
    result.imag_H_1 = result.imag_H_1 + sin(k*PI*(1-n)/50) + sin(k*PI*(5-n)/50);
    k++;
  }
  return result;
}

struct H_2 DFT_h_2(struct H_2 result, int n) // Discrete Fourier transform for h_2
{
  int k = 0;
  result.real_H_2 = 0;
  result.imag_H_2 = 0;
  while (k < 100)
  {
    result.real_H_2 = result.real_H_2 + h_2(k) * cos(k*PI*n/50);
    result.imag_H_2 = result.imag_H_2 + h_2(k) * -sin(k*PI*n/50);
    k++;
  }
  return result;
}

struct H_3 DFT_h_3(struct H_3 result, int n) // Discrete Fourier transform for h_3
{ 
  FILE *fp1; // declaring file handle for reading in h_3 data 
  fp1 = fopen("Values for h_3 against t.txt", "r");
  h_3 h3; 
  int k = 0;
  double t;
  result.real_H_3 = 0;
  result.imag_H_3 = 0;
  while (k < 200)
  {
    fscanf(fp1, "%i, %lf, %lf, %lf\n", &k, &t, &h3.real_h_3, &h3.imag_h_3);
    
    result.real_H_3 = result.real_H_3 + ((h3.real_h_3 * cos(k*PI*n/100)) + 
    (h3.imag_h_3 * sin(k*PI*n/100)));
      
    result.imag_H_3 = result.imag_H_3 + ((h3.imag_h_3 * cos(k*PI*n/100)) - 
    (h3.real_h_3 * sin(k*PI*n/100)));
    
    k++;
  }
  fclose(fp1);
  
  return result;
}

////////// Defining the values to be determined from applying inverse Fourier transforms ////////

typedef struct h_1_prime
{
  double real_h_1_prime, imag_h_1_prime;
} h_1_prime;

typedef struct h_2_prime
{
  double real_h_2_prime, imag_h_2_prime;
} h_2_prime;

typedef struct h_3_prime
{
  double real_h_3_prime, imag_h_3_prime;
} h_3_prime;

////////////////// Defining the inverse Fourier transforms using structs ////////////////////////

struct h_1_prime inv_DFT_H_1(struct h_1_prime result, int k)
{
  H_1 result_H1;
  int n = 0;
  result.real_h_1_prime = 0;
  result.imag_h_1_prime = 0;
  while(n < 100)
  { 
    if (n == 1) // omitting n = 1 from the summation by skipping to n = 2
    {
      n++;
    }
    else
    {
      result_H1 = DFT_h_1(result_H1, n);
      
      result.real_h_1_prime = result.real_h_1_prime + (result_H1.real_H_1 * cos(k*PI*n/50))/100;
      
      result.imag_h_1_prime = result.imag_h_1_prime + (result_H1.real_H_1 * sin(k*PI*n/50))/100; 
        // the second line here has real_H_1 because the imaginary part of the integral is still 
        //multiplied by a real number and the imaginary components of H_1 are all zero
      n++;
    }
  }
  return result;
}


struct h_2_prime inv_DFT_H_2(struct h_2_prime result, int k)
{
  H_2 result_H2;
  int n = 0;
  result.real_h_2_prime = 0;
  result.imag_h_2_prime = 0;
  while(n < 100)
  { 
    if (n == 0) // omitting n = 0 from the summation by skipping to n = 1
    {
      n++;
    }
    else
    {
      result_H2 = DFT_h_2(result_H2, n);
      
      result.real_h_2_prime = result.real_h_2_prime + (result_H2.real_H_2 * cos(k*PI*n/50))/100;
      
      result.imag_h_2_prime = result.imag_h_2_prime + (result_H2.real_H_2 * sin(k*PI*n/50))/100; 
        // only real part for the same reason mentioned in comment in inv DFT H_1
      n++;
    }
  }
  return result;
}

  
struct h_3_prime inv_DFT_H_3(struct h_3_prime result, int k)
{
  H_3 result_H3;
  int n = 0;
  result.real_h_3_prime = 0;
  result.imag_h_3_prime = 0;
  while(n < 200)
  {
    if (n == 1 || n == 3 || n == 4 || n == 13) // the four H_3 values with the largest magnitude
    {
    result_H3 = DFT_h_3(result_H3, n);
    
    result.real_h_3_prime = result.real_h_3_prime + ((result_H3.real_H_3 * cos(k*PI*n/100)) -
    (result_H3.imag_H_3 * sin(k*PI*n/100)))/200;
    
    result.imag_h_3_prime = result.imag_h_3_prime + ((result_H3.real_H_3 * sin(k*PI*n/100)) +
    (result_H3.imag_H_3 * cos(k*PI*n/100)))/200;
    
    n++;
    }
    else
    {
      n++;
    }
  }
  return result;
}


int main()
{
//////////////////////// Producing a text file for h_1 data /////////////////////////////////////

  int k; // defining loop variable
  
  FILE *fp2; // declare a file handle for h_1 data
  fp2 = fopen("Values for h_1 against t.txt", "w"); 

  for (k = 0; k < 100; k++)
  {
    fprintf(fp2, "%i, %.6e, %.6e, %.6e\n", k, k*PI/50, real_h_1(k), imag_h_1(k)); 
      // t values determined using expression involving k
  }
  fclose(fp2); //closes file for h_1 data
  
//////////////////////// Producing a text file for h_2 data /////////////////////////////////////

  FILE *fp3; // declares a file handle for h_2 data
  fp3 = fopen("Values for h_2 against t.txt", "w");
  
  for (k = 0; k < 100; k++) 
  {
    fprintf(fp3, "%i, %.6e, %.6e\n", k, k*PI/50, h_2(k)); 
  }
  fclose(fp3);

////////////////// Printing H_1 values using discrete Fourier transform /////////////////////////
  
  int n; // defining second loop variable
  
  struct H_1 *H1;
  H1 = (struct H_1*)malloc(100 * sizeof(struct H_1)); // assigning memory from the heap to store
    // store H_1 values
  
  for (n = 0; n < 100; n++)
  { 
    *H1 = DFT_h_1(*H1, n);
 
    printf("H_1 = %6.6f + %6.6fi for n = %i\n", H1->real_H_1, H1->imag_H_1, n);
  }
  free(H1); // deallocates memory used for H_1 values

////////////////// Printing H_2 values using discrete Fourier transform /////////////////////////

  struct H_2 *H2;
  H2 = (struct H_2*)malloc(100 * sizeof(struct H_2));

  for (n = 0; n < 100; n++)
  { 
    *H2 = DFT_h_2(*H2, n);
 
    printf("H_2 = %6.6f + %6.6fi for n = %i\n", H2->real_H_2, H2->imag_H_2, n);
  }
  free(H2);

//////////////////////// Producing a text file for h_1_prime values /////////////////////////////
  
  FILE *fp4; // declares a file handle for h_1_prime data
  fp4 = fopen("Values for h_1_prime against t.txt", "w");
  
  struct h_1_prime *h1prime;
  h1prime = (struct h_1_prime*)malloc(100 * sizeof(struct h_1_prime));
  
  for (k = 0; k < 100; k++)
  {
    *h1prime = inv_DFT_H_1(*h1prime, k);
    
    fprintf(fp4, " %i, %.6e, %.6e, %.6e\n", k, k*PI/50, h1prime->real_h_1_prime, 
    h1prime->imag_h_1_prime);
  }
  fclose(fp4);
  free(h1prime);
  
//////////////////////// Producing a text file for h_2_prime values /////////////////////////////

  FILE *fp5; // declares a file handle for h_2_prime data
  fp5 = fopen("Values for h_2_prime against t.txt", "w");
  
  struct h_2_prime *h2prime;
  h2prime = (struct h_2_prime*)malloc(100 * sizeof(struct h_2_prime));
  
  for (k = 0; k < 100; k++)
  {
    *h2prime = inv_DFT_H_2(*h2prime, k);
    
    fprintf(fp5, " %i, %.6e, %.6e, %.6e\n", k, k*PI/50, h2prime->real_h_2_prime, 
    h2prime->imag_h_2_prime);
  }
  fclose(fp5);
  free(h2prime);
  
//////////////////////// Producing a text file for h_3_prime values /////////////////////////////

  FILE *fp6; // declares a file handle for h_3_prime data
  fp6 = fopen("Values for h_3_prime against t.txt", "w"); 
  
  struct h_3_prime *h3prime;
  h3prime = (struct h_3_prime*)malloc(200 * sizeof(struct h_3_prime));
  
  for (k = 0; k < 200; k++)
  {  
    *h3prime = inv_DFT_H_3(*h3prime, k);
   
    fprintf(fp6, "%i, %.6e, %.6e, %.6e\n", k, k*PI/100, h3prime->real_h_3_prime, 
    h3prime->imag_h_3_prime);
  }
  fclose(fp6); 
  free(h3prime);
 
  return 0; // no errors to report
}

