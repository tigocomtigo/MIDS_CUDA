#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>
#define leitura 1
/*
__global__ void algoritmo(int* mat, int* res, int size, int* total_feromonios, int tamanho_menor_resposta){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int* marca = (int*) calloc(1,size * sizeof(int));
    int posicoes = size;

    int pos = rand() % size;
    int fer_ini = (rand() % total_feromonios[pos]) + 1;
    pos = -1;

    while(fer_ini > 0){
        ++pos;
        fer_ini -= mat[tid*size + pos];
    }

    marca[pos] = 1;
    --posicoes;

    for(int i = 0; i < size; ++i){
        if(mat[pos*size + i] == 0){
            marca[i] = 1;
            --posicoes;
        }
    }

    int controle;
    while(posicoes > 0){
        controle = 0;
        fer_ini = (rand() % total_feromonios[pos]) + 1;
        for(int i = 0; i < size; ++i){
            controle += mat[tid*size]
        }
    }

}
*/

__global__ void setup_kernel(curandState *state, unsigned long seed) {
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, id, 0, &state[id]);
}

__global__ void algoritmo(int* mat, int* res, int size, int* feromonios, int* menorGlobal, int* posicaoMenorGlobal, curandState *state, int* validacao, int* g_resFinal){
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int val = tid * size;
    if(tid < size){
        // Se aux == 0 a função deve fechar
        int aux = size;

        int menorLocal = 0;

        // Criando o vetor de validacao
        //int* validacao = (int*) malloc(size * sizeof(int));

        for(int i = 0; i < size; ++i){
            validacao[val + i] = 0;
            res[val + i] = 0;
        }

        // Seleciona aleatoriamente o primeiro vertice
        int pontape = curand(&state[tid]) % size;

        // Marcando que ja passou pelo pontape
        validacao[val + pontape] = 1;
        res[tid*size] = pontape;
        ++menorLocal;
        --aux;

        do{
            // Remover todos os valores invalidos
            for(int i = 0; i < size; ++i){
                if(mat[pontape*size + i] == 0 && validacao[val + i] != 1){
                    validacao[val + i] = 1;
                    --aux;
                }
            }

            // Seleciona aleatoriamente a primeira aresta
            int valor_caminhar_feromonios = curand(&state[tid]) % feromonios[pontape];

            int contador = 0;
            int valor_temporario = -1;
            for(int i = 0; i < size; ++i){
                contador += mat[pontape*size + i];
                if(validacao[val + i] == 0){
                    if(contador >= valor_caminhar_feromonios){
                        pontape = i;
                        validacao[val + pontape] = 1;
                        res[tid*size + menorLocal] = pontape;
                        ++menorLocal;
                        --aux;    
                        break;            
                    }else{
                        valor_temporario = i;
                    }
                }
                if(contador >= valor_caminhar_feromonios && valor_temporario != -1){
                    pontape = valor_temporario;
                    validacao[val + pontape] = 1;
                    res[tid*size + menorLocal] = pontape;
                    ++menorLocal;
                    --aux;
                    break;
                }
            }
        }while(aux > 0);

        res[tid*size + menorLocal] = -1;        

        if(menorLocal < menorGlobal[0]){
            menorGlobal[0] = menorLocal;
            posicaoMenorGlobal[0] = tid;
            int i = 0;
            while(res[tid*size + i] != -1){
                g_resFinal[i] = res[tid*size + i];
                ++i;
            }
            g_resFinal[i] = -1;
            //res[tid*size + menorLocal] = -1;
        }

    }
    //__syncthreads();
}

__global__ void atualizaFeromonios(int* mat, int* res, int* posMenor, int size, int* feromonios){
    int pivo = 0;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid < size){
        while(res[pivo + 1] != -1){
            mat[res[pivo]*size + res[pivo + 1]] += 20;
            mat[res[pivo + 1]*size + res[pivo]] += 20;
            feromonios[res[pivo]] += 20;
            feromonios[res[pivo + 1]] += 20;
            ++pivo;
        }        
    }

}

int main(){


    FILE *file;
    file = fopen("instancia_teste.mis", "r");
    int size, vertices;

    srand(time(NULL));
    clock_t start, end;
    double cpu_time_used;

    if (leitura == 1)
    {
        fscanf(file, "%i %i\n", &size, &size);
        fscanf(file, "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ");
    }
    else if (leitura == 2)
    {
        fscanf(file, "p edge %i %i    ", &size, &vertices);
    }

    int* c_mat = (int*) malloc(size * size * sizeof(int));
    int* c_res = (int*) malloc(size * size * sizeof(int));
    int* feromoniosVetor = (int*) calloc(1, size * sizeof(int));
    for(int i = 0; i < size; ++i){
        //total_feromonios = 0;
        for(int j = 0; j < size; ++j){
            c_mat[(i*size) + j] = 1;
            //++total_feromonios;
        }
        //feromoniosVetor[i] = total_feromonios;
    }

    int valor1, valor2;

    if (leitura == 1){
        while (!feof(file))
        {
            fscanf(file, "%i %i\n", &valor1, &valor2);
            // printf("%i %i\n", valor1, valor2);
            //  fscanf(file,"%i %i\n", &valor1, &valor2);
            --valor1;
            --valor2;
            c_mat[valor1*size + valor2] = 0;
            c_mat[valor2*size + valor1] = 0;
            ++feromoniosVetor[valor1];
            ++feromoniosVetor[valor2];
        }
    }else if(leitura == 2){

    }

    int* menor = (int*) malloc(sizeof(int));
    menor[0] = size;
    int* posMenor = (int*) malloc(sizeof(int));
    posMenor[0] = 0;

    int* g_mat;
    int* g_res;
    int* g_fer;
    int* g_menor;
    int* g_posMenor;
    int* g_validacao;

    cudaMalloc(&g_mat, size * size * sizeof(int));
    cudaMalloc(&g_res, size * size * sizeof(int));
    cudaMalloc(&g_fer, size * sizeof(int));
    cudaMalloc(&g_menor, sizeof(int));
    cudaMalloc(&g_posMenor, sizeof(int));
    cudaMalloc(&g_validacao, size * size * sizeof(int));

    cudaMemcpy(g_mat, c_mat, size * size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_fer, feromoniosVetor, size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_menor, menor, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_posMenor, posMenor, sizeof(int), cudaMemcpyHostToDevice);

    //int blockSize = 512;
    //int numBlocks = (size + blockSize - 1)/blockSize;

    // Numeros aleatorios
    curandState *dev_state;

    // Aloca memória na GPU
    cudaMalloc((void **)&dev_state, size * sizeof(curandState));

    // Configura o estado inicial do gerador de números aleatórios
    setup_kernel<<<1, size>>>(dev_state, time(NULL));

    int* g_res_final;
    cudaMalloc(&g_res_final, size * sizeof(int));

    int ciclos = 100;
    start = clock();
    while(ciclos > 0){
        algoritmo<<<1, size>>>(g_mat, g_res, size, g_fer, g_menor, g_posMenor, dev_state, g_validacao, g_res_final);
        atualizaFeromonios<<<1, size>>>(g_mat, g_res_final, g_posMenor, size, g_fer);
        --ciclos;
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    int* c_res_final = (int*) malloc(size * sizeof(int));

    cudaMemcpy(menor, g_menor, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(posMenor, g_posMenor, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_res_final, g_res_final, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_res, g_res,size * size * sizeof(int), cudaMemcpyDeviceToHost);

    printf("%i\n", menor[0]);
    printf("%i\n", posMenor[0]);
    /*
    for(int i = 0; i < size; ++i){
        printf(" %i", c_res[posMenor[0]*size + i]);
    }    
    */
   for(int i = 0 ; i < size; ++i){
        if(c_res_final[i] == -1){
            break;
        }
        printf("%i ", c_res_final[i] + 1);
        //printf("\n");
   }
    printf("\nOutros resultados:\n");

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if(c_res[i*size + j] == -1){
                break;
            }
            printf("%i ", c_res[i*size + j] + 1);
        }
        printf("\n");
    }

    printf("Tempo de execução: %.8f segundos\n", cpu_time_used);

    /*
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            printf("%i ", c_mat[i*size + j]);
        }
        printf("\n");
    }    
    */

   cudaFree(g_fer);
   cudaFree(g_mat);
   cudaFree(g_menor);
   cudaFree(g_posMenor);
   cudaFree(g_res);
   cudaFree(g_validacao);

   free(c_mat);
   free(c_res);
   free(feromoniosVetor);
   free(menor);
   free(posMenor);

    return 0;
}