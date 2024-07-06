/*-------------------------------------------------------------------------
 *              UNIFAL - UNIVERSIDADE Federal de Alfenas.
 *                  BACHEARELADO EM CIENCIA DA COMPUTACAO.
 * Trabalho..:  Contagem de feijoes - Complementar
 * Professor.:  Luiz Eduardo da Silva
 * Aluno.....:  Pedro Henrique Alves Barbosa - 2022.1.08.043
 * Data......:  20/05/2024
 -------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RED_PIXEL_VALUE 255

#if defined(_WIN32) || defined(__WIN64__) || defined(__CYGWIN__)
#include "..\\utils\\imagelib.h"
#elif defined(__linux__)
#include "../utils/imagelib.h"
#endif

typedef struct no *ptno;
typedef struct no {
    int i, j;
    ptno next;
} no;

/*Insere um novo nó na fila circular. Se a fila está vazia, o novo nó aponta para si mesmo. 
Caso contrário, o novo nó é inserido logo após o último nó da fila.*/
void insQ(ptno *Q, int i, int j) {
    ptno new = malloc(sizeof(no));
    new->i = i;
    new->j = j;
    if (!(*Q))
        new->next = new;
    else {
        new->next = (*Q)->next;
        (*Q)->next = new;
    }
    *Q = new;
}

/*Remove o nó da frente da fila circular. A função ajusta os ponteiros para manter a integridade da fila 
e atualiza os valores i e j com os dados do nó removido.*/
void remQ(ptno *Q, int *i, int *j) {
    ptno init = (*Q)->next;
    *i = init->i;
    *j = init->j;
    if (*Q == init)
        *Q = NULL;
    else
        (*Q)->next = init->next;
    free(init);
}

/*Verifica se a fila está vazia. Retorna 1 se a fila estiver vazia, caso contrário, retorna 0.*/
int isEmpty(ptno *Q) {
    return *Q == NULL;
}
/*Inicializa uma fila de prioridade com mn níveis. Cada nível é uma fila circular que pode conter nós.*/
void initQPrior(ptno *QPrior, int mn) {
    int i;
    for (i = 0; i < mn; i++)
        QPrior[i] = NULL;
}

/*Insere um elemento na fila de prioridade, usando a função insQ para adicionar o elemento na fila do nível p.*/
void insert(ptno *QPrior, int i, int j, int p) {
    insQ(QPrior + p, i, j);
}

/*Remove o próximo elemento da fila de prioridade. A função busca o elemento de maior prioridade (menor índice) 
disponível e atualiza os valores i e j com os dados do nó removido. Se todas as filas estiverem vazias, retorna 1; 
caso contrário, retorna 0.*/
int pop(ptno *QPrior, int *i, int *j, int *maxPrior, int mn) {
    while (*maxPrior < mn && isEmpty(QPrior + *maxPrior))
        (*maxPrior)++;
    if (*maxPrior == mn)
        return 1;
    remQ(QPrior + *maxPrior, i, j);
    return 0;
}

/*Calcula o gradiente de intensidade da imagem In usando uma janela de raio raio. O valor do gradiente em cada pixel 
é a diferença entre o valor máximo e mínimo dos pixels dentro da janela definida pelo raio.*/
image gradient(image In, int raio) {
    int i, j, y, x, max, min;
    int nl = In->nr, nc = In->nc, mn = In->ml;
    image Out = img_clone(In);
    for (i = 0; i < nl; i++)
        for (j = 0; j < nc; j++) {
            max = -1;
            min = mn + 1;
            for (y = -raio; y <= raio; y++)
                for (x = -raio; x <= raio; x++) {
                    int pi = i + y;
                    int pj = j + x;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc) {
                        if (abs(x) + abs(y) <= raio && In->px[pi * nc + pj] > max)
                            max = In->px[pi * nc + pj];
                        if (abs(x) + abs(y) <= raio && In->px[pi * nc + pj] < min)
                            min = In->px[pi * nc + pj];
                    }
                }
            Out->px[i * nc + j] = max - min;
        }
    return Out;
}

/*Aplica o algoritmo de watershed na imagem In, começando a partir da posição (y, x). 
O algoritmo segmenta a imagem identificando regiões de interesse (contornos), que são 
marcadas em vermelho na imagem de saída Out.*/
image watershed(image In, int y, int x) {
    int i, j, k, maxPrior = 0, stop = 0;
    int nl = In->nr, nc = In->nc, mn = In->ml;
    ptno qPrior[mn];
    image mark = img_create(nl, nc, mn, In->tp);
    image Out = img_create(nl, nc, mn, In->tp);

    enum {
        NONE, QUEUE, WSHED, MARK1, MARK2,
    };

    struct {
        int i, j;
    } n4[4] = { {0, 1}, {1, 0}, {0, -1}, {-1, 0} };

    initQPrior(qPrior, mn);

    for (i = 0; i < nl * nc; i++)
        mark->px[i] = NONE;

    int raio = 10;
    for (i = -raio; i <= raio; i++)
        for (j = -raio; j <= raio; j++) {
            int pi = i + y;
            int pj = j + x;
            if (abs(i) + abs(j) <= raio)
                mark->px[pi * nc + pj] = MARK1;
        }

    for (i = 0; i < nl; i++) {
        mark->px[i * nc] = MARK2;
        mark->px[i * nc + nc - 1] = MARK2;
    }
    for (j = 0; j < nc; j++) {
        mark->px[j] = MARK2;
        mark->px[(nl - 1) * nc + j] = MARK2;
    }

    for (i = 1; i < nl - 1; i++)
        for (j = 1; j < nc - 1; j++)
            if (mark->px[i * nc + j] == NONE) {
                int isAdj = 0;
                for (k = 0; k < 4; k++) {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    int m = mark->px[pi * nc + pj];
                    if (m == MARK1 || m == MARK2)
                        isAdj = 1;
                }
                if (isAdj) {
                    mark->px[i * nc + j] = QUEUE;
                    insert(qPrior, i, j, In->px[i * nc + j]);
                }
            }

    while (!stop) {
        int m = NONE;
        int isWshed = 0;
        stop = pop(qPrior, &i, &j, &maxPrior, mn);
        if (!stop) {
            for (k = 0; k < 4; k++) {
                int pi = i + n4[k].i;
                int pj = j + n4[k].j;
                if (pi >= 0 && pi < nl && pj >= 0 && pj < nc) {
                    int mAdj = mark->px[pi * nc + pj];
                    if (mAdj == MARK1 || mAdj == MARK2) {
                        if (m == NONE)
                            m = mAdj;
                        else if (m != mAdj)
                            isWshed = 1;
                    }
                }
            }
            if (isWshed)
                mark->px[i * nc + j] = WSHED;
            else {
                mark->px[i * nc + j] = m;
                for (k = 0; k < 4; k++) {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc)
                        if (mark->px[pi * nc + pj] == NONE) {
                            int prior, px;
                            mark->px[pi * nc + pj] = QUEUE;
                            px = In->px[pi * nc + pj];
                            prior = (px < maxPrior) ? maxPrior : px;
                            insert(qPrior, pi, pj, prior);
                        }
                }
            }
        }
    }
    // Marcar os contornos vermelhos na imagem de saída
    for (i = 0; i < nl; i++) {
      for (j = 0; j < nc; j++) {
        if (mark->px[i * nc + j] == WSHED && In->px[i * nc + j] > 160) {
            Out->px[i * nc + j] = RED_PIXEL_VALUE;  // Marcar como contorno vermelho
        } else {
            Out->px[i * nc + j] = In->px[i * nc + j];  // Copiar pixels não modificados
        }
    }
}
    img_free(mark);
    return Out;
}


// Função para calcular o negativo de uma imagem
image neg_pgm(image In) {
    image Out = img_clone(In);
    for (int i = 0; i < In->nr * In->nc; i++)
        Out->px[i] = In->ml - In->px[i];
    return Out;
}

// Função para aplicar um limiar à imagem de entrada
image intensidade(image In) {
    float Limiar[In->ml + 1];
    image Out = img_clone(In);
    for (int i = 0; i < In->ml + 1; i++) {
        Limiar[i] = i < 160 ? 0 : 1;
    }
    for (int i = 0; i < In->nr * In->nc; i++)
        Out->px[i] = Limiar[In->px[i]];
    return Out;
}

// Função similar à intensidade, mas com outro limiar
image intensidade2(image In) {
    float Limiar[In->ml + 1];
    image Out = img_clone(In);
    for (int i = 0; i < In->ml + 1; i++) {
        Limiar[i] = i < 6 ? 0 : 1;
    }
    for (int i = 0; i < In->nr * In->nc; i++)
        Out->px[i] = Limiar[In->px[i]];
    return Out;
}

int find(int parent[], int i) {
    while (parent[i] != i)
        i = parent[i];
    return i;
}

void Union(int parent[], int i, int j) {
    int x = find(parent, i);
    int y = find(parent, j);
    parent[y] = x;
}

// Função para contar o número de rótulos únicos em uma imagem rotulada
int contar_rotulos(image img, int rotulos[]) {
    int cont = 0;
    int rotuloUsado[img->nc * img->nr];
    for (int i = 0; i < img->nc * img->nr; i++)
        rotuloUsado[i] = 0;
    for (int i = 0; i < img->nc * img->nr; i++) {
        int rotuloRaiz = find(rotulos, img->px[i]);
        if (rotuloUsado[rotuloRaiz] == 0 && rotuloRaiz != 0) {
            rotuloUsado[rotuloRaiz] = 1;
            cont++;
        }
    }
    return cont;
}

// Função para rotular os componentes conectados em uma imagem binária
void label(image In) {
    int numeroLinhas = In->nr;
    int numeroColunas = In->nc;
    int *pixels = In->px;
    int numRotulo = 0;
    int parent[1000];

    // Inicializa o array parent com índices
    for (int i = 0; i < 1000; i++)
        parent[i] = i;

    // Percorre a imagem para rotular os componentes
    for (int i = 1; i < numeroLinhas; i++) {
        for (int j = 1; j < numeroColunas; j++) {
            int p = pixels[i * numeroColunas + j];
            int r = pixels[(i - 1) * numeroColunas + j];
            int t = pixels[i * numeroColunas + j - 1];
            if (p != 0) {
                if (r == 0 && t == 0)
                    pixels[i * numeroColunas + j] = ++numRotulo;
                if (r != 0 && t == 0)
                    pixels[i * numeroColunas + j] = r;
                if (r == 0 && t != 0)
                    pixels[i * numeroColunas + j] = t;
                if (r != 0 && t != 0 && t == r)
                    pixels[i * numeroColunas + j] = r;
                if (r != 0 && t != 0 && t != r) {
                    pixels[i * numeroColunas + j] = t;
                    Union(parent, r, t);
                }
            }
        }
    }

    // Atualiza os rótulos dos pixels com os seus pais
    for (int i = 0; i < numeroLinhas * numeroColunas; i++)
        In->px[i] = find(parent, In->px[i]);

    // Define a quantidade total de componentes conectados
    In->ml = numRotulo;

    // Imprime o número de componentes conectados
    printf("#componentes= %d\n", contar_rotulos(In, parent));
}

// Função para encontrar o mínimo entre três números
int minimo3(int a, int b, int c) {
    if (a < b && a < c) {
        return a;
    }
    if (b < c) {
        return b;
    }
    return c;
}

// Função para calcular a distância de cada pixel ao fundo da imagem binária
int distancia(image In) {
    int maxDist = -1;
    int nLinhas = In->nr;
    int nColunas = In->nc;
    int *pixels = In->px;

    // Calcula a distância da esquerda para a direita e de cima para baixo
    for (int i = 1; i < nLinhas - 1; i++) {
        for (int j = 1; j < nColunas - 1; j++) {
            int p = pixels[i * nColunas + j];
            int a = pixels[(i - 1) * nColunas + j];
            int b = pixels[i * nColunas + j - 1];

            if (p != 0) {
                pixels[i * nColunas + j] = (a + 1) < (b + 1) ? (a + 1) : (b + 1);
            }
        }
    }

    // Calcula a distância da direita para a esquerda e de baixo para cima
    for (int i = nLinhas - 2; i > 0; i--) {
        for (int j = nColunas - 2; j > 0; j--) {
            int p = pixels[i * nColunas + j];
            int a = pixels[i * nColunas + j + 1];
            int b = pixels[(i + 1) * nColunas + j];

            if (p != 0) {
                pixels[i * nColunas + j] = minimo3(a + 1, b + 1, p);
                if (maxDist < pixels[i * nColunas + j]) {
                    maxDist = pixels[i * nColunas + j];
                }
            }
        }
    }
    return maxDist;
}

// Função para exibir uma mensagem de uso do programa
void exibir_mensagem(char *s) {
    printf("\nContafeijao");
    printf("\n-------------------------------");
    printf("\nUso:  %s  nome-imagem[.pgm] \n\n", s);
    printf("    nome-imagem[.pgm] é o nome do arquivo da imagem \n");
    exit(1);
}

int main(int argc, char *argv[]) {
    char nameIn[100], nameOut[100], comando[110];
    image In, Out, Grd;

    // Verifica se o número de argumentos é válido
    if (argc < 2)
        exibir_mensagem(argv[0]);

    // Obtém os nomes dos arquivos de In e Out
    img_name(argv[1], nameIn, nameOut, GRAY, GRAY);  // Alterando para COLOR para imagem de saída

    // Lê a imagem de In
    In = img_get(nameIn, GRAY);

    // Aplica as transformações na imagem
    Out = neg_pgm(In);
    Out = intensidade(Out);
    Out->ml = distancia(Out);
    Out = intensidade2(Out);
    label(Out);

    // Aplica a função watershed
    Grd = gradient(Out, 1);
    Out = watershed(Grd, In->nr / 2, In->nc / 2);

    // Salva a imagem resultante
    img_put(Out, nameOut, COLOR);

    // Exibe a imagem resultante
    sprintf(comando, "%s %s &", VIEW, nameOut);
    system(comando);

    // Libera a memória alocada
    img_free(In);
    img_free(Out);

    return 0;
}