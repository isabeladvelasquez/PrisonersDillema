/***************************************************************************************
*                        Dilema do Prisioneiro - Densidade Variável                    *
*                                         24/11/2025                                   *
***************************************************************************************/
//gcc -o 'L100_RHOvar_M0' 'L100_RHOvar_M0.c' -I/home/isabela/IC/lat2eps/lat2eps-master/ -O3 -lgsl -lgslcblas -lm


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "lat2eps.h"

/***************************************************************************************
*                                Random Numbers Functions                              *
*                                                                                      *
***************************************************************************************/

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

/****************************************************************************************
*                                  Definição de Constantes                              *
*                                                                                       *
****************************************************************************************/

#define L 100                 //sites in each line of the net
#define N (L * L)             //sites in the net
#define MCS 10000             //maximum simulation time (MCS)

#define S 0.0                 //pontuação de C quando encontra D
#define TT 1.4                //pontuação de D quando encontra C
#define R 1.0                 //pontuação de C quando encontra C
#define P 0.0                 //pontuação de D quando encontra D

#define SAVE_CONFIG 0         //salvar a configuração final
#define SAVE_TEMP 0           //salvar pc em relação ao tempo
#define SAVE_ATVSITE 0        //save the last config of the net showing active sites

/***************************************************************************************
*                               Declaring other Functions                              *
*                                                                                      *
***************************************************************************************/

//Rede do sistema
float rho; //densidade teórica da rede
float n; //número de agentes na rede
float p[21]; //densidade real total da rede
float pcp[21]; //densidade de cooperadores pela densidade da rede
float pcp_tempo[21][44]; //densidade de cooperadores pela densidade da rede em relação ao tempo

//Time in MCS
int t[21][44];
int timetarget[] = {1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20, 24, 29, 35, 43, 52, 63, 77, 94, 115, 140, 170, 207, 252, 307, 374, 455, 554, 674, 820, 998, 1214, 1477, 1798, 2187, 2661, 3237, 3938, 4791, 5829, 7091, 8627, 10000};
int savet = 0;

//Points each agent has on the round
float pontos[N];

//Estratégias e posições iniciais
int s[N]; //s[i] = 0 se for D, s[i] = 1 se for C e s[i] = 2 se o sítio for vazio
int vizinhanca[N][4]; /*para cada i em N temos 4 vizinhos onde: [0] = direita; [1] = esquerda; [2] = acima; [3] = abaixo.*/
int s_atv[N];

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
void pts_iniciais(){for (int i = 0; i < N; ++i){pontos[i] = 0.0;}} //todos começam com 0 pontos em cada rodada

//Estratégias em função da densidade da rede
void estrategias(){
    n = 0.0;
    for (int i = 0; i < N; i++){
        if (randi() <= rho) {
          n += 1;
          if (randi() < 0.5) s[i] = 0;
          else s[i] = 1;} 
        else {s[i] = 2; } } }

//Densidade de C em N
float densidade(){
    float soma_C = 0.0;
    for (int k = 0; k < N; k++) {
        if (s[k] == 1) {soma_C += 1.0;}
    }
    return soma_C/N;}

/******************************************************************************************
*                               FUNÇÕES DE ARQUIVOS DE DADOS                              *
*                                                                                         *
******************************************************************************************/

//Função para a seed inicial
unsigned long funcSeed(){
  unsigned long seed = (unsigned long)time(NULL);
  if (seed%2 == 0) seed += 1; //semente é ímpar
  return seed;
  }

FILE*fp1;
char name[100];
int ok = 0;
//Salvar os dados de pc e pcp
void saveConfig(unsigned long seed){
  if (SAVE_CONFIG == 1){
    while (ok == 0){
      sprintf(name,"SDP_L%d_RHO(0-1)_M0_S%lu.dat",L,seed);
      FILE *temp_fp = fopen(name, "r");
      if (temp_fp != NULL){seed += 1; fclose(temp_fp); ok = 0;}
      else ok = 1;}
    fp1 = fopen(name,"w");
    fprintf(fp1,"# Simulação do Dilema do Prisioneiro\n");
    fprintf(fp1,"# Comprimento da rede: %d\n",L);
    fprintf(fp1,"# RHO: (0-1)\n");
    fprintf(fp1,"# Tempo máximo de simulação: %d\n", MCS);
    fprintf(fp1,"# seed rho p pcp\n"); }
  else fp1 = stdout;}

//Salvar pcp_tempo e tempo
int saveTemp(unsigned long seed){
  if (SAVE_TEMP == 1){
    while (ok == 0){
      sprintf(name,"SDP_L%d_RHO(0-1)_M0_S%lu.dat",L,seed);
      FILE *temp_fp = fopen(name, "r");
      if (temp_fp != NULL){seed += 1; fclose(temp_fp); ok = 0;}
      else ok = 1;}
    fp1 = fopen(name,"w");
    fprintf(fp1,"# Simulação do Dilema do Prisioneiro - rho em relação ao tempo\n");
    fprintf(fp1,"# Comprimento da rede: %d\n",L);
    fprintf(fp1,"# RHO: %.2f\n", rho);
    fprintf(fp1,"# Tempo máximo de simulação: %d\n", MCS);
    fprintf(fp1,"# seed rho p t pcp\n"); }
  else fp1 = stdout;}

//Net graph
void plot_atvsite(int time){
if (SAVE_ATVSITE == 1){
    snprintf(name, sizeof(name), "SDP(ATVS)_L100_RHO%.2f_M0_t%d.eps", rho, time);
    lat2eps_init(L,L);
    lat2eps_set_color(1,0x0115b2); // blue
    lat2eps_set_color(0,0xFF0000); // red
    lat2eps_set_color(2,0xFFFFFF); // white
    lat2eps_set_color(3,0xdbc400); // yellow
    
    for (int j=0; j<N; ++j) {
        if (s[j] == 1) lat2eps_set_site(j%L,j/L,1);
        if (s[j] == 0) lat2eps_set_site(j%L,j/L,0);
        if (s[j] == 2) lat2eps_set_site(j%L,j/L,2);
        if (s_atv[j] == 1) lat2eps_set_site(j%L,j/L,3);
        }

    lat2eps_gen_eps(name,0,0,L,L,2,6);

    lat2eps_release();
}}

/*****************************************************************************************************************************
*                                                        Simulação                                                           *
*****************************************************************************************************************************/
int main(){
  unsigned long seed = funcSeed(); //define a semente de números aleatórios
    
  saveConfig(seed); //Funções de registro de dados
  saveTemp(seed);
  
  init_rng(seed);
  
  
  p[0] = 0;
  pcp[0] = 0;
  for (int i = 0; i < 44; ++i){t[0][i] = 0; pcp_tempo[0][i] = 0;}
  rede(); //montando a rede

/********************************* Interações ******************************************************************************/
  for (int l = 1; l < 21; ++l){ //loop para determinar a densidade teórica
    rho = l/20.0;
    estrategias();
    savet = 0;
    for (int i = 0; i < N; ++i) {s_atv[i] = 0;}
    plot_atvsite(0);
    for (int k = 0; k < MCS; ++k) { //loop para as interações
      pts_iniciais();
      for (int i = 0; i < N; ++i){ //jogo de cada jogador N
      
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
            if (s[vizinhanca[i][j]] == 1){ //neighbor agente is C
              pontos[i] += R;}
          }
        }
      }
      
    int s_novo[N]; //atualização de estratégias ao final de cada rodada vvvv
    for (int i = 0; i < N; ++i){
      if (s[i] < 2){
        int agente = i;
        int j0 = randi() * 4;
        float ponto = pontos[i];
        for (int j = 0; j < 4; ++j){
          if (ponto < pontos[vizinhanca[i][j0]]) {
            ponto = pontos[vizinhanca[i][j0]];
            agente = vizinhanca[i][j0];}
          j0 = (j0 + 1) % 4;}
        s_novo[i] = s[agente]; }
      else s_novo[i] = 2; }
      
    for (int i = 0; i < N; ++i) {
      if (s_novo[i] != s[i]) s_atv[i] = 1; else s_atv[i] = 0;
      s[i] = s_novo[i];} 
    
    p[l] = n/N;
    if(SAVE_TEMP == 1 && (k + 1) == timetarget[savet]){
      t[l][savet] = k + 1;
      pcp_tempo[l][savet] = densidade()/p[l];
      savet += 1; } } 
    
/*********************************** Recolhendo dados no final da simulação **************************************************/
  plot_atvsite(10000);
  if (SAVE_CONFIG == 1){pcp[l] = densidade()/p[l];  fprintf(fp1,"%lu %.4f %.4f %.4f\n", seed, rho, p[l], pcp[l]);}
  if (SAVE_TEMP == 1){for (int j = 0; j < 44; ++j) fprintf(fp1,"%lu %.4f %.4f %d %.4f\n", seed, rho, p[l], t[l][j], pcp_tempo[l][j]);}
  }

  if (SAVE_CONFIG == 1 || SAVE_TEMP == 1) fclose(fp1);
  liberar_rng();
  return 0;
}
