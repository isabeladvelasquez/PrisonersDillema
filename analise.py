import sys
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import expit

########################################################################################################################
########################################## Reading the Archive's Lines #################################################
########################################################################################################################

def auxi(fileName, column):
    linhas_limpas = []
    try:
        with open(fileName, 'r') as f:
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            dados_linear = f.read().split()
            
        for i in range(0, len(dados_linear), column):
            try:
                bloco_str = dados_linear[i : i + column]
                if len(bloco_str) == column:
                    bloco_float = [float(d) for d in bloco_str]
                    linhas_limpas.append(bloco_float)
            except (ValueError, IndexError):
                continue
    except FileNotFoundError:
        return None
    except Exception as e:
        return None
        
    if not linhas_limpas:
        return None
    return np.array(linhas_limpas)

########################################################################################################################
################################################## Calculating Means ###################################################
########################################################################################################################

def process(lines):
    path = os.path.join(folder, file)
    listaArquivo = glob.glob(path)
    
    listaDados = []
    for arquivo in listaArquivo:
        dadosArquivo = auxi(arquivo, column)
        
        if dadosArquivo is not None and dadosArquivo.shape[0] == lines:
            listaDados.append(dadosArquivo)
    
    if not listaDados:
        return None, None
        
    dados_3d = np.stack(listaDados, axis=0)
    
    mediaTotal = np.mean(dados_3d, axis=0)
    media_x = mediaTotal[:, column_1_index] 
    media_y = mediaTotal[:, column_2_index] 
    return media_x, media_y, dados_3d

########################################################################################################################
################################################# PLOT FUNCTIONS #######################################################
########################################################################################################################

def verifica_cprob(dados):
    # --- INÍCIO DA VERIFICAÇÃO ADICIONADA ---
    # Coluna de índice 4 em todas as simulações e todas as linhas
    coluna_4 = dados[:, :, 4] 
    
    # Verifica se os valores são 0.0 ou 1.0 (com margem de erro mínima)
    mask_erro = ~np.logical_or(np.isclose(coluna_4, 0.0, atol=1e-7), 
                                np.isclose(coluna_4, 1.0, atol=1e-7))

    if np.any(mask_erro):
        valores_estranhos = np.unique(coluna_4[mask_erro])
        print(f"Foram encontrados {len(valores_estranhos)} valores decimais na coluna 4!")
        print(f"Valores que não são 0 ou 1: {valores_estranhos}")
    return

def plot_prob_c(x, y, dados, samples, n_rhos, n_m):
    colors = plt.cm.turbo(np.linspace(0, 1, n_m))
    c_prob = y/samples
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1,2,1)
    for i in range(n_m):
        cor = colors[i]
        m = (i+1)/n_m
        start = i * n_rhos
        end = start + n_rhos

        plt.plot(x[start:end], y[start:end], color = cor, label = r'$m' f'= {m:.1f}$')
    
    plt.title(r'Prob($\rho_c/\rho = 1$) vs $\rho$ with mobility')
    plt.xlim(0, 1)
    plt.ylabel(r'Prob($\rho_c/\rho = 1$)')
    plt.xlabel(r'$\rho$')
    plt.ylim(0, 0.2)
    plt.tight_layout()
    plt.legend(loc='upper left', fontsize=10, title=r"Mobility $m$") 
    
    plt.subplot(1,2,2)
    for j in [0, 4, 9]:
        cor = colors[j]
        m = (j+1)/n_m
        start = j * n_rhos
        end = start + n_rhos

        plt.plot(x[start:end], y[start:end], color = cor, label = r'$m' f'= {m:.1f}$')
    
    plt.title(r'Prob($\rho_c/\rho = 1$) vs $\rho$ with mobility (m = 0.1, 0.5, 1.0)')
    plt.xlim(0, 1)
    plt.ylabel(r'Prob($\rho_c/\rho = 1$)')
    plt.xlabel(r'$\rho$')
    plt.ylim(0, 0.2)
    plt.tight_layout()
    plt.legend(loc='upper left', fontsize=10, title=r"Mobility $m$") 
    plt.savefig(f"PROB2C_VS_RHO.pdf")
    plt.show()


def plot_TIMEvsPCP_MGRAFS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #pra 10 rhos
    ajuste = 1 #pra todos os rhos
    cores = plt.cm.turbo(np.linspace(0, 1, int(n_rhos/ajuste)))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    for i in range(n_m):
        m = (i+1)/n_m
        plt.figure(1, figsize=(10, 10))
        for j in range(int(n_rhos/ajuste)): 
            rho = (j+1)/(n_rhos/ajuste)
            inicio = (i * n_m * n_time) + ajuste*j * n_time
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]
            
            for k in range(samples):
                plt.plot(dados_organizados[k, i, ajuste*j, :, column_1_index], dados_organizados[k, i, ajuste*j, :, column_2_index], ls = '-', color=cores[j], alpha = 0.3, linewidth = 0.8)

            plt.plot(xconfig, yconfig, ls = '-', color=cores[j], label = r'$\rho$' f' = {rho}') 

        plt.xscale('log')
        plt.title(f"Densidade de Cooperadores entre os Agentes pelo Tempo com Mobilidade m = {m}")
        plt.xlim(1, 10000)
        #plt.ylim(0, 0.7)
        plt.xlabel(r"$t$(MCS)")
        plt.ylabel(r"$\rho_{c}/\rho$") 
        
        plt.legend(ncol=2, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title=r"Densidade $\rho$")
        plt.tight_layout()
        
        plt.savefig(f"graf_TIMEVSPCP_M{i+1}_10RHOS.pdf" if ajuste != 1 else f"graf_TIMEVSPCP_M{i+1}.pdf")
        plt.close()


def plot_TIMEvsPCP_RHOGRAFS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #pra 10 rhos
    ajuste = 1 #pra todos os rhos
    cores = plt.cm.turbo(np.linspace(0, 1, int(n_m)))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    for i in range(int(n_rhos/ajuste)):
        rho = (i+1)/(n_rhos/ajuste)
        plt.figure(1, figsize=(10, 10))
        for j in range(n_m): 
            m = (j+1)/n_m
            inicio = (ajuste-1)*n_time + (i * n_time) + (j * n_rhos * n_time)
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]
            
            for k in range(samples):
                plt.plot(dados_organizados[k, j, ajuste*i, :, column_1_index],dados_organizados[k, j, ajuste*i, :, column_2_index], ls = '-', color=cores[j], alpha = 0.3, linewidth = 0.8)

            plt.plot(xconfig, yconfig, ls = '-', color=cores[j], label = r'$m$' f' = {m}') 

        plt.xscale('log')
        plt.title(f"Densidade de Cooperadores entre os Agentes pelo Tempo com Densidade = {rho}")
        plt.xlim(1, 10000)
        #plt.ylim(0, 0.7)
        plt.xlabel(r"$t$(MCS)")
        plt.ylabel(r"$\rho_{c}/\rho$") 
        
        plt.legend(ncol=2, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title=r"Mobilidade $m$")
        plt.tight_layout() # Ajusta as margens

        plt.savefig(f"graf_TIMEVSPCP_RHO{ajuste*(i+1)}.pdf" if ajuste != 1 else f"graf_TIMEVSPCP_RHO{i+1}.pdf")
        plt.close()


def plot_TIMEvsPCP_RHOGRAFS_INDIVIDUAL(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    cores = plt.cm.turbo(np.linspace(0, 1, n_m))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    for i in range(n_rhos): 
        rho = (i+1)/n_rhos
        for j in range(n_m):
            m = (j+1)/n_m
            plt.figure(figsize=(10, 7))
            
            for k in range(samples):
                plt.plot(dados_organizados[k, j, i, :, column_1_index], dados_organizados[k, j, i, :, column_2_index], ls = '-', color=cores[j], alpha = 0.15, linewidth = 0.8)
                
            inicio = (i * n_time) + (j * n_rhos * n_time)
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]
            
            plt.plot(xconfig, yconfig, ls = '-', color=cores[j], label = r'$m' f' = {m:.2f}$')
            plt.xscale('log')
            plt.title(r'Densidade de Cooperadores pelo Tempo para $m$ =' f'{(j+1)}'r' e $\rho$ =' f'{rho:.3f}')
            plt.xlim(1, 10000)
            #plt.ylim(0, 0.7)
            plt.xlabel(r"$t$(MCS)")
            plt.ylabel(r"$\rho_C / \rho$") 
            
            plt.legend()
            plt.tight_layout() # Ajusta as margens
            
            plt.savefig(f"grafico_TIMEVSPCP_INDIVIDUAL_M{j+1}_RHO{i+1}.pdf")
            plt.close()


def plot_TIMEvsPCP_RHOGRAFS_MEANS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #pra 10 rhos
    ajuste = 1 #pra todos os rhos
    cores = plt.cm.turbo(np.linspace(0, 1, n_m))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    for i in range(int(n_rhos/ajuste)):
        rho = (i+1)/(n_rhos/ajuste)
        plt.figure(1, figsize=(10, 7))
        for j in range(n_m): 
            m = (j+1)/n_m
            inicio = (ajuste-1)*n_time + (i * n_time) + (j * n_rhos * n_time)
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]

            plt.plot(xconfig, yconfig, ls = '-', color=cores[j], label = r'$m$' f' = {m:.1f}') 

        plt.xscale('log')
        plt.title(f"Densidade de Cooperadores entre os Agentes pelo Tempo com Densidade = {rho:.3f}")
        plt.xlim(1, 10000)
        #plt.ylim(0, 0.7)
        plt.xlabel(r"$t$(MCS)")
        plt.ylabel(r"$\rho_{c}/\rho$") 
        
        plt.legend(ncol=2, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title=r"Mobilidade $m$")
        plt.tight_layout() 
        
        plt.savefig(f"graf_TIMEVSPCP_RHO{i+1}_MEAN_10RHOS.pdf" if ajuste != 1 else f"graf_TIMEVSPCP_RHO{i+1}_MEAN.pdf")
        plt.close()


def plot_TIMEvsPCP_MGRAFS_MEANS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #pra 10 rhos
    ajuste = 1 #pra todos os rhos
    cores = plt.cm.turbo(np.linspace(0, 1, int(n_rhos/ajuste)))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    for i in range(n_m):
        m = (i+1)/n_m
        plt.figure(1, figsize=(10, 10))
        for j in range(int(n_rhos/ajuste)): 
            rho = (j+1)/(n_rhos/ajuste)
            inicio = (i * n_m * n_time) + ajuste*j * n_time
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]

            plt.plot(xconfig, yconfig, ls = '-', color=cores[j], label = r'$\rho$' f' = {rho}') 

        plt.xscale('log')
        plt.title(f"Densidade de Cooperadores entre os Agentes pelo Tempo com Mobilidade m = {m}")
        plt.xlim(1, 10000)
        #plt.ylim(0, 0.7)
        plt.xlabel(r"$t$(MCS)")
        plt.ylabel(r"$\rho_{c}/\rho$") 
        
        plt.legend(ncol=2, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title=r"Densidade $\rho$")
        plt.tight_layout() # Ajusta as margens
        
        plt.savefig(f"graf_TIMEVSPCP_M{i+1}_MEAN_10RHOS.pdf" if ajuste != 1 else f"graf_TIMEVSPCP_M{i+1}_MEAN.pdf")
        plt.close()


def plot_3d_TIMEVSPCPVSM(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #para 10 rhos
    ajuste = 1
    cores = plt.cm.turbo(np.linspace(0, 1, int(n_rhos/ajuste)))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    fig = plt.figure(figsize=(13, 11))
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(n_m):
        m = (i+1)/n_m
        plt.figure(1, figsize=(10, 10))
        for j in range(int(n_rhos/ajuste)): 
            rho = (j+1)/(n_rhos/ajuste)
            inicio = (i * n_m * n_time) + ajuste*j * n_time
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]
            
            legenda = r'$\rho = $'f'{rho}' if i == 0 else None
            ax.plot(xconfig, yconfig, zs=i, zdir='y', c=cores[j], label=legenda)

    # Make legend, set axes limits and labels
    ax.legend(ncol=4, bbox_to_anchor=(1.3, 0.8), loc='right', borderpad=1)
    ax.set_xlim(0, 10000)
    ax.set_ylim(0, 10)
    ax.set_zlim(0, 0.2)
    ax.set_xlabel('T (MCS)')
    ax.set_ylabel('M')
    ax.set_zlabel(r'$\rho_C / \rho$')
    ax.grid(False)
    ax.set_box_aspect((3, 4.8, 3))

    ax.view_init(elev=20., azim=-140, roll=0)
    plt.savefig('3d_TIMEVSPCPVSM_10RHOS.png' if ajuste != 1 else f"3d_TIMEVSPCPVSM.pdf", dpi=300, bbox_inches='tight')
    plt.show()


def plot_3d_TIMEVSPCPVSP(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time):
    #ajuste = int(n_rhos/10) #para 10 rhos
    ajuste = 1 #para todos rhos
    cores = plt.cm.turbo(np.linspace(0, 1, n_m))
    dados_organizados = dados.reshape(samples, n_m, n_rhos, n_time, column)
    
    fig = plt.figure(figsize=(13, 11))
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(int(n_rhos/ajuste)):
        rho = (i+1)/(n_rhos/ajuste)
        plt.figure(1, figsize=(10, 10))
        for j in range(n_m): 
            m = (j+1)/n_m
            inicio = (ajuste-1)*n_time + (i * n_time) + (j * n_rhos * n_time)
            fim = inicio + n_time
            
            xconfig = x[inicio:fim]
            yconfig = y[inicio:fim]
            
            legenda = r'$m = $'f' {m}' if i == 0 else None
            ax.plot(xconfig, yconfig, zs=i+1, zdir='y', c=cores[j], label=legenda)

    # Make legend, set axes limits and labels
    ax.legend(ncol=4, bbox_to_anchor=(1.3, 0.8), loc='right', borderpad=1)
    ax.set_xlim(0, 10000)
    ax.set_ylim(0, int(n_rhos/ajuste))
    ax.set_zlim(0, 0.2)
    ax.set_xlabel('T (MCS)')
    ax.set_ylabel(r'$\rho$')
    ax.set_zlabel(r'$\rho_C / \rho$')
    ax.grid(False)
    ax.set_box_aspect((3, 4.8, 3))

    ax.view_init(elev=20., azim=-140, roll=0)
    plt.savefig('3d_TIMEVSPCPVSP_10RHOS.png' if ajuste != 1 else f"3d_TIMEVSPCPVSP.pdf", dpi=300, bbox_inches='tight')
    plt.show()

########################################################################################################################
##################################################### PLOTS ############################################################ 
########################################################################################################################

file = "*.dat"
folder = "/home/isabela/IC/2CPROB_TIMEVSPCP_062026/SPD_2C_L100_T10000_2CPROB_062026"
samples = 100

column = 6
lines = 400

column_1_index = 3
column_2_index = 5

n_time = 76
n_rhos = 40
n_m = 10

if __name__ == "__main__":
    x, y, dados = process(lines)

    if x is not None and y is not None and dados is not None:

        try:
            plt.rcParams.update({
                "text.usetex": True,
                "font.family": "serif",
                "font.serif": ["Computer Modern Roman"],
                "font.size": 14
            })
        except RuntimeError:
            pass 
        #FUNCTIONS THAT WILL BE USED BELOW
        
        #plot_prob_c(x, y, dados, samples, n_rhos, n_m)
        #verifica_cprob(dados)

        ############## TIMEVSPCP WITH SAMPLES ############
        #plot_TIMEvsPCP_MGRAFS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)  
        #plot_TIMEvsPCP_RHOGRAFS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)   
        #plot_TIMEvsPCP_RHOGRAFS_INDIVIDUAL(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)

        ############## TIMEVSPCP WITHOUT SAMPLES #########
        #plot_TIMEvsPCP_RHOGRAFS_MEANS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)
        #plot_TIMEvsPCP_MGRAFS_MEANS(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)
        #plot_3d_TIMEVSPCPVSM(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)
        #plot_3d_TIMEVSPCPVSP(x, y, dados, samples, column, column_1_index, column_2_index, n_rhos, n_m, n_time)
    else:
        pass
        
    
