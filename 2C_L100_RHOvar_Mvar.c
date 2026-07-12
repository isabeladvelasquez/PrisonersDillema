/***************************************************************************************************************
*                            Prisoner's Dillema Simulation - Variable Density and Motion 2C                    *
*                                                   28/06/2026                                                 *
***************************************************************************************************************/
//gcc '2C_L100_RHOvar_Mvar.c' lat2eps_lib.c -o '2C_L100_RHOvar_Mvar' -O3 -lgsl -lgslcblas -lm

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "lat2eps.h"

/*******************************************************************************************************************
*                                                RANDOM NUMBERS GENERATOR                                          *
*                                                                                                                  *
*******************************************************************************************************************/

//GSL
gsl_rng *r;
void init_rng(unsigned long int seed) {
    const gsl_rng_type * T = gsl_rng_mt19937; 
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
}
void liberar_rng() {gsl_rng_free(r);}

//Random number between 0 and 1
float randi() {return (float)gsl_rng_uniform(r);}

/**********************************************************************************************************************
*                                                   DEFINING CONSTANTS                                                *
*                                                                                                                     *
**********************************************************************************************************************/
#define L 100                 //sites in each line of the net
#define N (L * L)             //sites in the net
#define MCS 10000             //maximum simulation time (MCS)

#define S 0.0                 //pontuação de C quando encontra D (sucker)
#define TT 1.4                //pontuação de D quando encontra C (temptation)
#define R 1.0                 //pontuação de C quando encontra C (reward)
#define P 0.0                 //pontuação de D quando encontra D (punishment)

#define SAVE_CONFIG 0         //salvar a configuração final
#define SAVE_TIME 1           //salvar pc em relação ao tempo
#define SAVE_ATVSITE 1        //save the last config of the net showing active sites
#define SAVE_2CPROB 0         //save the probability of invasion by cooperators

#define NRHOS 40              //number of rho stamps
#define NTIMETARGETS 76       //number of time stamps

/***************************************************************************************************************************
*                                                  DECLARING FUNCTIONS                                                     *
*                                                                                                                          *
***************************************************************************************************************************/

//System's net
float rho; //theoretical network density (variable)
double rhos[NRHOS]; //theoretical network density (array)
double create_rhos(){ //generate an array of rhos
  for (int n = 0; n < NRHOS; ++n){rhos[n] = (double)(n+1)/NRHOS;} }
float n; //number of agents on the network
float p[NRHOS]; //real network density
float pcp[NRHOS]; //cooperators density by the density of the network
float pcp_tempo[NRHOS][NTIMETARGETS]; //cooperators density by the density of the network during time
int c_prob; //probability of invasion by C (it is only determined at the end of the simulation)

//Time in MCS
int timetarget[NTIMETARGETS];
void create_timetargets(){
  double exponent = log10(MCS)/NTIMETARGETS;
  int last_value = 0;
  for (int n = 0; n < NTIMETARGETS; ++n){
    int value = round(pow(10, (double)(n+1) * exponent));
    while (value <= last_value){value += 1;}
    last_value = value;
    timetarget[n] = last_value; } }
int savet = 0;

//Motion and movement
float m;

//Points each agent has on the round
float pontos[N];

//Initial strategies and positions
int s[N]; //s[i] = 0 se for D, s[i] = 1 se for C e s[i] = 2 se o sítio for vazio
int vizinhanca[N][4]; /*para cada i em N temos 4 vizinhos onde: [0] = direita; [1] = esquerda; [2] = acima; [3] = abaixo.*/
int s_atv[N];
int s_novo[N];
int end_sim; //if == 0, the system has achieved a stationary state

//Neighborhood net
void rede(){
    for (int i = 0; i < N; ++i) {
            if (i % L == L - 1) 
                vizinhanca[i][0] = i + 1 - L;
            else 
                vizinhanca[i][0] = i + 1;

            if (i % L == 0) 
                vizinhanca[i][1] = i - 1 + L;
            else 
                vizinhanca[i][1] = i - 1;

            if (i < L) 
                vizinhanca[i][2] = i - L + N;
            else 
                vizinhanca[i][2] = i - L;

            if (i >= N - L) 
                vizinhanca[i][3] = (i % L);
            else 
                vizinhanca[i][3] = i + L;
        }
    }

//Initial points
void init_points(){for (int i = 0; i < N; ++i) pontos[i] = 0.0;} //every agent starts the round with 0 points

//Strategies based on the density of the network
void Cstrategy(){
    n = 0.0;
    for (int i = 0; i < N; i++) {
        if (randi() <= rho) {
          n += 1;
          s[i] = 0;}
        else {s[i] = 2;} }
    int x = gsl_rng_uniform_int(r, N);
    if (s[x] == 2) {n += 1;}
    s[x] = 1;
    if (s[vizinhanca[x][0]] == 2) {n += 1;}
    s[vizinhanca[x][0]] = 1; }

//Densidade de C em N
float densidade(){
    float soma_C = 0.0;
    for (int k = 0; k < N; k++) {
        if (s[k] == 1) {soma_C += 1.0;}
    }
    return soma_C/N;}

//Setting all active sites to zero
void zero_atvsites(){memset(s_atv, 0, N * sizeof(int)); }

/**********************************************************************************************************************
*                                                   DATA SAVING FUNTIONS                                              *
*                                                                                                                     *
**********************************************************************************************************************/

//Função para a seed inicial
unsigned long funcSeed(){
  unsigned long seed = (unsigned long)time(NULL);
  if (seed%2 == 0) seed += 1; //semente é ímpar
  return seed;
  }

FILE *fp1, *fp2, *fp3;
char name1[100]; char name2[100]; char name3[100];
int ok;
//Salvar os dados de pc e pcp
void save(unsigned long seed){
  if (SAVE_CONFIG == 1){
    ok = 0;
    while (ok == 0){
      sprintf(name1,"SPD_2C_L%d_T%d_SAVECONFIG_S%lu.dat",L,MCS,seed);
      FILE *temp_fp1 = fopen(name1, "r");
      if (temp_fp1 != NULL){seed += 1; fclose(temp_fp1); ok = 0;}
      else ok = 1;}
    fp1 = fopen(name1,"w");
    fprintf(fp1,"# Prisoner's Dillema Simulation - rho_c vs rho\n");
    fprintf(fp1,"# Comprimento da rede: %d\n",L);
    fprintf(fp1,"# Saving system's final 'rho_c/rho' configuration after simulation\n");
    fprintf(fp1,"# Maximum simulation time: %d\n", MCS);
    fprintf(fp1,"# seed m rho p pcp\n");} 
  else fp1 = stdout;
  if (SAVE_TIME == 1) {
    ok = 0;
    while (ok == 0){
      sprintf(name2,"SPD_2C_L%d_T%d_SAVETIME_S%lu.dat",L,MCS,seed);
      FILE *temp_fp2 = fopen(name2, "r");
      if (temp_fp2 != NULL){seed += 1; fclose(temp_fp2); ok = 0;}
      else ok = 1;}
    fp2 = fopen(name2,"w");
    fprintf(fp2,"# Prisoner's Dillema Simulation - rho vs time\n");
    fprintf(fp2,"# Network length: %d\n",L);
    fprintf(fp2,"# Saving system's 'rho_c/rho' configuration during in simulation steps\n");
    fprintf(fp2,"# Maximum simulation time: %d (MCS)\n", MCS);
    fprintf(fp2,"# seed m rho p t pcp\n"); }
  else fp2 = stdout;
  if (SAVE_2CPROB == 1) {
    ok = 0;
    while (ok == 0){
      sprintf(name3,"SPD_2C_L%d_T%d_SAVE2CPROB_S%lu.dat",L,MCS,seed);
      FILE *temp_fp3 = fopen(name3, "r");
      if (temp_fp3 != NULL){seed += 1; fclose(temp_fp3); ok = 0;}
      else ok = 1;}
    fp3 = fopen(name3, "w");
    fprintf(fp3,"# Prisoner's Dillema Simulation - probability of C invasion vs rho\n");
    fprintf(fp3,"# Network's length: %d\n",L);
    fprintf(fp3,"# Saving system's variable 'c_prob' if C invades at the end of the simulation\n");
    fprintf(fp3,"# Maximum simulation time: %d (MCS)\n", MCS);
    fprintf(fp3,"# seed m rho p pcp cprob\n"); }
  else fp3 = stdout; }

char name[256];
//Net graph
void plot_atvsite(float time, unsigned long seed){
if (SAVE_ATVSITE == 1){
    snprintf(name, sizeof(name), "SDP_ATVS_L100_RHO%.3f_M%.1f_T%.5f_S%lu.eps", rho, m, time, seed);
    lat2eps_init(L,L);
    lat2eps_set_color(1,0x0115b2); // blue
    lat2eps_set_color(0,0xFF0000); // red
    lat2eps_set_color(2,0xFFFFFF); // white
    lat2eps_set_color(3,0xdbc400); // yellow
    lat2eps_set_color(4,0xff38c3); // pink
    lat2eps_set_color(5,0x05a815); // green
    
    for (int j=0; j<N; ++j) {
        if (s[j] == 1) lat2eps_set_site(j%L,j/L,1);
        if (s[j] == 0) lat2eps_set_site(j%L,j/L,0);
        if (s[j] == 2) lat2eps_set_site(j%L,j/L,2);
        if (s_atv[j] == 1 && s[j] == 1) lat2eps_set_site(j%L,j/L,3);
        if (s_atv[j] == 1 && s[j] == 0) lat2eps_set_site(j%L,j/L,4);
        }

    lat2eps_gen_eps(name,0,0,L,L,2,6);

    lat2eps_release();
}}

/****************************************************************************************************************************
*                                                                                                                           *
*                                                        SIMULATION                                                         *
*                                                                                                                           *
****************************************************************************************************************************/

int main(){
  unsigned long seed = funcSeed(); //define a semente de números aleatórios
    
  save(seed); //Funções de registro de dados
  
  init_rng(seed);
  
  rede(); 
  create_rhos();
  create_timetargets();
  /*************************************************** STEP 1: COMBATING ***************************************************/
  for (int count_m = 6; count_m < 7 ; ++count_m){ //movement loop
    m = count_m/10.0;
    for (int l = 15; l < 16; ++l){ //theoretical density loop
      rho = rhos[l];
      Cstrategy();
      plot_atvsite(0, seed); 
      c_prob = 0;
      zero_atvsites();
      end_sim = 0;
      savet = 0;
      for (int k = 0; k < MCS; ++k) { //rounds loop
        init_points(); 
        zero_atvsites();
        if (end_sim == 0) {
          end_sim = 1;
          for (int i = 0; i < N; ++i){ //turn of each agent N
          
            if (s[i] == 0) { //agent chosen is D
              for (int j = 0; j < 4; ++j){
                if (s[vizinhanca[i][j]] == 0){ //neighbor agent is D
                  pontos[i] += P;}
                if (s[vizinhanca[i][j]] == 1){ //neighbor agent is C
                      pontos[i] += TT;}
              }
            }
            if (s[i] == 1){ //agent chosen is C
              for (int j = 0; j < 4; ++j){
                if (s[vizinhanca[i][j]] == 0){ //neighbor agent is D
                  pontos[i] += S;}
                if (s[vizinhanca[i][j]] == 1){ //neighbor agent is C
                  pontos[i] += R;}
              }
            }
          }
      
  /**************************************************** STEP 2: DIFUSION **************************************************/
          if (m != 0.0){          
            for (int d = 0; d < N; ++d){
              int i = gsl_rng_uniform_int(r, N);
              if (s[i] != 2) {
                int j0 = gsl_rng_uniform_int(r, 4);
                if (s[vizinhanca[i][j0]] == 2 && randi() <= m){
                  s[vizinhanca[i][j0]] = s[i]; 
                  pontos[vizinhanca[i][j0]] = pontos[i]; 
                  s[i] = 2;
                  pontos[i] = 0; 
                  if ((s[vizinhanca[i][j0]] == 1) && (k < 10)) plot_atvsite(((float)k) + (((float)(i + 1)) / (2*N)), seed); } } } 
          } 

  /***************************************** STEP 3: OFFSPRING AND ACTIVE SITES *******************************************/  
          memcpy(s_novo, s, N * sizeof(int));    
   

          for (int i = 0; i < N; ++i){
            if (s[i] != 2){
              int agente = i;
              int j0 = gsl_rng_uniform_int(r, 4);
              float ponto = pontos[i];
              for (int j = 0; j < 4; ++j){
                if (ponto < pontos[vizinhanca[i][j0]]) {
                  ponto = pontos[vizinhanca[i][j0]];
                  agente = vizinhanca[i][j0];}
                j0 = (j0 + 1) % 4;} 
              s_novo[i] = s[agente]; } }

          if (SAVE_ATVSITE == 1) {
            for (int i = 0; i < N; ++i) {
              if (s_novo[i] != s[i]) s_atv[i] = 1; 
              else s_atv[i] = 0;
          } }

          memcpy(s, s_novo, N * sizeof(int));
    
/********************************************* SAVING DATA FOR A SINGLE ROUND **********************************************/
          p[l] = n/N;
          pcp[l] = densidade()/p[l];
          if(SAVE_TIME == 1 && savet < NTIMETARGETS && (k + 1) == timetarget[savet]){
            pcp_tempo[l][savet] = pcp[l];
            savet += 1; } 
          if (pcp[l] == 0 || pcp[l] == 1) end_sim = 1;
          else end_sim = 0; }

        else {
          if(SAVE_TIME == 1 && savet < NTIMETARGETS && (k + 1) == timetarget[savet]){
            pcp_tempo[l][savet] = pcp[l];
            savet += 1; } } } 
        
      plot_atvsite(MCS, seed); 
      if (pcp[l] == 1.0)c_prob = 1;
/******************************************* SAVING DATA FOR A SINGLE RHO AND M ********************************************/
      
      if (SAVE_CONFIG == 1){fprintf(fp1,"%lu %.1f %.4f %.4f %.4f\n", seed, m, rho, p[l], pcp[l]);}
      if (SAVE_TIME == 1){for (int j = 0; j < NTIMETARGETS; ++j) fprintf(fp2,"%lu %.1f %.4f %.4f %d %.4f\n", seed, m, rho, p[l], timetarget[j], pcp_tempo[l][j]); }
      if (SAVE_2CPROB == 1){fprintf(fp3,"%lu %.1f %.4f %.4f %.4f %d\n", seed, m, rho, p[l], pcp[l], c_prob);} 
      fflush(fp1); fflush(fp2); fflush(fp3);} }

/************************************************** ENDING SIMULATION ******************************************************/
  if (SAVE_CONFIG == 1) fclose(fp1);
  if (SAVE_TIME == 1) fclose(fp2);
  if (SAVE_2CPROB == 1) fclose(fp3);
  liberar_rng();
  return 0;
  
}
