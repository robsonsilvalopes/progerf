
/******************************************************************************

data    : 19 de Novembro de 2013
Horário : 14:10
Autor   : Robson da Silva Lopes
Alteracao: gerar degenerações

GPAVLIMCFYWHKRQNEDST
0	G 
1	P 
2	A 
3	V 
4	L 
5	I 
6	M 
7	C 
8	F 
9	Y 
10	W 
11	H 
12	K 
13	R 
14	Q 
15	N 
16	E 
17	D 
18	S 
19	T 

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LAST(k,n) ((k) & ((1<<(n))-1)) // pega os n bits menos significativos de k inteiro

typedef struct lista_fragmentos
{
    char * frag; //Guarda os caracteres do bloco repetitivo
    int  lenght; //Guarda a quantidade de repetições
    int posicaoini; //Guarda a posicao inicial do bloco repetitivo
    int posicaofim; //Guarda a posicao final do bloco repetitivo
    int gaps;	// gaps existentes no bloco;
    char * deg;
    //int quant; //Guarda a quantidade fragmento anteriores a ele dentro de um bloco repetitivo
    struct lista_fragmentos * pro; //Aponta para o proximo fragmento
    struct lista_fragmentos * ant; //Aponta para o fragmento anterior
} lista;

typedef struct lista_degen{
    char        * frag; //Guarda os caracteres do bloco repetitivo
    unsigned long int         hashvalue;      // armazena o valor hash
    struct lista_degen * proximo; // aponta para o proximo elemento.
    struct lista_degen * anterior; //aposta para o elemento anteior
} listaDegen;

typedef struct no_cabeca
{
    lista * first; //Aponta para o primeiro bloco repetitivo em uma posicao da tabela hash
    lista * last; //Aponta para o ultimo fragmento em uma posicao da tabela hash
} cabeca;

typedef struct no_head
{
    listaDegen * first;
    int contVariation;
    listaDegen * last;
    
} head;

lista           * aloca         (char * elem, int posini, int posfim, int janela);
listaDegen      * alocaDegen    (char * elem, int valorHash, int janela);
unsigned long int             funhash         (char * frag, int janela);
__inline unsigned long int    funhashSimplified(char a);
unsigned long int             insere          (cabeca ** hash, char * elem, int pos, int janela, int rep, int gap, long inserir);
unsigned long int             insereDegen     (head ** hashDegen,char * elem, int janela, long inserir);
void            insereDegenVariation(head ** hashDegen,char * elem, int janela, int position);
void            insereDegeneration(cabeca ** hash, head ** hashDegen, char * elem, int pos, int janela, int rep, int gap, long inserir);
void            generate        (head ** hashDegen,char * elem, int iniPos, int janela, int txDegen, int deep, int position);
void            libera          (cabeca ** hash, unsigned long int tam);
int             nulo            (char * elem, int janela);
__inline int    nuloChar(char * elem, int janela);
void            resultado       (FILE * fp, cabeca ** hash, char * nome, int tam, int janela, int rep);
void            liberaHashDegen (head ** hash, unsigned long int tam);
fpos_t          pos; //armazena a posição do stream que varre o arquivo fpin
char            nucleotideos[] = "GPAVLIMCFYWHKRQNEDST";

int main(int argc, char ** argv)
{
    FILE * fpin, *fpout;
    int count, f, i, janela, rep, gap,flag, flagEntryDegen,qtNuclChange;
    int countAux;
    unsigned long int p;
    long antHash = 0;
    float txDegeneration =  20;
    char c, *elem, *in, *nome, *out;
    cabeca ** hash;
    head ** hashDegen; // Função hash que armazena as degenerações possíveis
    if (argc < 9) {
        puts("Modo de execucao invalido.");
        puts("A chamada do programa de ser da seguinte maneira:");
        puts("<nome_programa> -q <arquivo_fasta> -o <arquivo_de_saida> -i <tamanho_da_subsequencia> -t <tamanho_minimo_do_bloco_repetitivo>");
        puts("ex <nome_programa> -q sequence.fasta -o saida.ssr -i 4 -t 5 -g 0 -d 1");
        exit(0); //Sai do programa
    }
    for (i = 1; i < argc; i = i + 2) {
        if (!strcmp("-q", argv[i])) {
            in          = (char*) calloc(strlen(argv[i + 1]), sizeof (char));
             strcpy(in, argv[i + 1]); //Pega o nome do arquivo de entrada
         } else if (!strcmp("-o", argv[i])) {
            out = (char*) calloc(strlen(argv[i + 1]), sizeof (char));
            strcpy(out, argv[i + 1]); //Pega o nome do arquivo de saida
        } else if (!strcmp("-i", argv[i])) {
            janela = atoi(argv[i + 1]); //Pega o tamanho que cada fragmento tera
        } else if (!strcmp("-t", argv[i])) {
            rep = atoi(argv[i + 1]); //Pega o numero minimo de repetiticoes necessirias
        } else if (!strcmp("-g", argv[i])) {
            gap = atoi(argv[i + 1]); //Pega o numero minimo de gaps
        } else if (!strcmp("-d", argv[i])) {
            txDegeneration = atoi(argv[i + 1]); //Pega a taxa de degeneração
        }
        else {
            fprintf(stderr, "Paramento invilido.\n");
            exit(1); //Sai do programa
        }
    }
    if (janela < 2 || janela > 12) {
        fprintf(stderr, "Tamanho da subsequencia invalido.\n");
        exit(4);
    }
    if ((fpin = fopen(in, "r")) == NULL) { //Tenta abrir o arquivo de entrada
        fprintf(stderr, "Erro ao abrir arquivo de leitura.\n");
        exit(2); //Sai do programa
    }
    if ((fpout = fopen(out, "w")) == NULL) { //Tenta abrir o arquivo de saida
        fprintf(stderr, "Erro ao abrir arquivo de escrita.\n");
        exit(3); //Sai do programa
    }
    
    
    p           = pow(20, janela); //Calcula o tamanho da tabela hash
    hash        = (cabeca**) calloc(p, sizeof (cabeca*)); //Aloca a tabela hash
    hashDegen   = (head**) calloc(p, sizeof (head*)); //Aloca a tabela hash
    elem        = (char*) calloc(janela, sizeof (char)); //Aloca o fragmento
    nome        = (char*) calloc(30, sizeof (char)); //Aloca o vetor que conteri o ID da sequencia genetica
    qtNuclChange = janela*(txDegeneration/100);
    
    //realiza o pré-processamento gerando a hash com as degenerações possíveis
    f           = 2; //Nao seri necessiria imprimir dados da tabela hash, ou seja, primeira sequencia lida
    flag        = 1; // flag que define se vai armazenar a posição da stream
    flagEntryDegen = 0;
    while ((fscanf(fpin, "%c", &c)) != EOF) {
        if (c == '>') { //Inicio da sequencia genetica
            
            flag                = 1;                         
            flagEntryDegen      = 1;
            
            if (f != 2) { //Nao e a primeira sequencia lida do arquivo de entrada
                resultado(fpout, hash, nome, count + janela - 1, janela, rep); //Imprime resultado para o arquivo de saida
                libera(hash, p); //Apaga os dados da tabela hash
                liberaHashDegen(hashDegen,p); // apaga os dados da tabela hash

            }
            f = 1; //Na proxima interacao seri lido o primeiro caracter da sequencia
            fscanf(fpin, "%s", nome); //Pega o ID da sequencia genetica
            while (c != '\n') fscanf(fpin, "%c", &c); //Le o restante da linha que contem o ID
            
            
            if(flag){
                fgetpos (fpin,&pos); //Retrieves the current position in the stream.
                flag = 0;
            }
        } else if (c != '\n') { //'Enters' no meio da sequencia
            
            
            
            if(flagEntryDegen){
                
                if(f){
                    f = 0 ;
                    count = 0; //Inicia o contador do tamanho da sequencia
                    elem[0] = toupper(c); //Guada o primeiro caracter na primeira posicao fragmento convertendo-o para maiusculo
                    for (i = 1; i < janela; i++) { //Leitura do restantes de carecteres para preenche o fragmento
                        fscanf(fpin, "%c", &c); //Leitura do um caracter
                        elem[i] = toupper(c); //Converte para maiusculo e guarda na proxima posicao do fragmento
                    }
                    if (nulo(elem, janela)) //Verifica se nao contem nucleotideos nao identificados
                    {
 //                       printf("\t%i \t- %s\n", count, elem);
                        antHash = insereDegen(hashDegen, elem, janela, -1);
                    }
                }else{
                    memmove(&elem[0],&elem[1],(janela -1)* sizeof(char));
                    elem[janela - 1] = toupper(c);
                    if(nulo(elem, janela))
                    {
 //                       printf("\t%i \t- %s\n", count, elem);
                        antHash = insereDegen(hashDegen, elem, janela, antHash);
                    }
                }
            }
             
             
                         
            if (f && !flagEntryDegen) { //Leitura do primeira caracter da sequencia genetica
                f = 0; //Para leitura de apenas um caracter
                count = 0; //Inicia o contador do tamanho da sequencia
                elem[0] = toupper(c); //Guada o primeiro caracter na primeira posicao fragmento convertendo-o para maiusculo
                for (i = 1; i < janela; i++) { //Leitura do restantes de carecteres para preenche o fragmento
                    fscanf(fpin, "%c", &c); //Leitura do um caracter
                    elem[i] = toupper(c); //Converte para maiusculo e guarda na proxima posicao do fragmento
                }
                
                if (nulo(elem, janela)) //Verifica se nao contem proteinas nao identificadas
                    antHash = insere(hash, elem, count, janela, rep, gap, -1); //Insercao na tabela hash
            } else if(!f && !flagEntryDegen ) {
                memmove(&elem[0], &elem[1], (janela - 1) * sizeof (char)); //Desloca os (janela - 1) caracteres uma posicao para a esquerda
                elem[janela - 1] = toupper(c); //Converte para maiusculo e guarda na ultima posicao do fragmento
                if (nuloChar(elem, janela)){ //Verifica se nao contem proteinas nao identificadas
                    insereDegeneration(hash,hashDegen,elem,count,janela,rep,gap,antHash); // verifica a degeneração e insere no hash de soluçao
                    antHash = insere(hash, elem, count, janela, rep, gap, antHash); //Insercao na tabela hash
                    
                    
                } 
                   
            }
            count++; //Incrementa o contador de caracteres.
        }else{
            if(flagEntryDegen) // inserir combinações possíveis;
            {
                for(countAux = 0; countAux < p; countAux++ ) // percorre todas a combinações possíveis
                {
                   
                    if(hashDegen[countAux] != NULL){
                         generate(hashDegen,hashDegen[countAux]->first->frag,0,janela,qtNuclChange,0, countAux);
                    }
                }
                
                fsetpos (fpin,&pos);
                flagEntryDegen = 0;
                f = 1;
            }
            
        }
    }
    resultado(fpout, hash, nome, count + janela - 1, janela, rep); //Imprime o resultado da ulima sequencia do arquivo de entrada
    libera(hash, p); //Apaga os dados da tabela hash
    fclose(fpin); //Fecha o arquivo de entrada
    fclose(fpout); //Fecha o arquivo de saida
    free(hash); //Libera memoria da tabela hash
    free(hashDegen); //Libera memoria da tabela hashDegen
}

lista * aloca(char * elem, int posini, int posfim, int janela)
{
 // printf("aloca \n");
    lista * novo = (lista*) calloc(1, sizeof (lista)); //Cria uma nova estruta
    novo->frag = (char*) calloc(janela, sizeof (char)); //Cria um vetor dentro da estrutura
    strcpy(novo->frag, elem); //Copia fragmento para dentro da estrura
    novo->posicaoini = posini; //Seta posicao inicial do fragmento
    novo->posicaofim = posfim;
    novo->gaps = 0;
    novo->lenght = 1;
    novo->deg = (char*) calloc(2, sizeof (char)); //Cria um vetor dentro da estrutura
    novo->deg = "";
    
 //   printf("out aloca\n");
    return novo;
}


listaDegen * alocaDegen(char * elem, int valorHash, int janela)
{
 // printf("alocaDegen \n");
    listaDegen * novo = (listaDegen*) calloc(1, sizeof (listaDegen)); //cria uma nova estrutura
    novo->frag     = (char*) calloc(janela, sizeof (char)); //Cria um vetor dentro da estrutura
    strcpy(novo->frag, elem);
    novo->hashvalue = valorHash;
 //   printf("out alocaDegen\n");
    return novo;
    
}
unsigned long int funhash(char * frag, int janela)
{
  // cada letra será representado por um número na base 19
  // para repetiticoes em DNA utilizou-se a base binária
   // printf("funhash \n");
    int i, t;
    unsigned long int ret = 0;
    t = janela - 1;
    for (i = 0; i < janela; i++) {
        switch (frag[i]) {
	  // A letra G = 0
        case 'P':	// P 1
            ret += (pow(19, (t - i)));
            break;
        case 'A':	// A 2
            ret += (2*pow(19, (t - i)));
            break;
	case 'V':	// V 3
            ret += (3*pow(19, (t - i)) );
            break;
	case 'L':	// L 4
            ret += (4*pow(19, (t - i)));
            break;
	case 'I':	// I 5
            ret += (5*pow(19, (t - i)));
            break;
	case 'M':	// M 6
            ret += (6*pow(19, (t - i)));
            break;
	case 'C':	// C 7
            ret += (7*pow(19, (t - i)));
            break;
	case 'F':	// F 8
            ret += (8*pow(19, (t - i)));
            break;
	case 'Y': 	// Y 9
            ret += (9*pow(19, (t - i)));
            break;
	case 'W':	// W 10
            ret += (10*pow(19, (t - i)));
            break;
	case 'H':	// H 11
            ret += (11*pow(19, (t - i)));
            break;
	case 'K':	// K 12
            ret += (12*pow(19, (t - i)));
            break;
	case 'R':	// R 13
            ret += (13*pow(19, (t - i)));
            break;
        case 'Q':	// Q 14
            ret += (14*pow(19, (t - i)));
            break;
	case 'N':	// N 15
            ret += (15*pow(19, (t - i)));
            break;
	case 'E':	// E 16
            ret += (16*pow(19, (t - i)));
            break;
	case 'D':	// D 17
            ret += (17*pow(19, (t - i)));
            break;
	case 'S':	// S 18
            ret += (18*pow(19, (t - i)));
            break;
        case 'T':	// T 19
            ret += (19*pow(19, (t - i)));
	
        }
    }
//    printf("out funhash\n");
    return ret;
}

__inline unsigned long int funhashSimplified(char a){
  
    
    int ret = 0;
//    printf("funhashSimplified \n");
    switch(a){
      // A letra G = 0
        case 'P':
            ret = 1;
            break;
        case 'A':
            ret = 2;
            break;
	case 'V':
            ret = 3;
            break;
	case 'L':
            ret = 4;
            break;
	case 'I':
            ret = 5;
            break;
	case 'M':
            ret = 6;
            break;
	case 'C':
            ret = 7;
            break;
	case 'F':
            ret = 8;
            break;
	case 'Y':
            ret = 9;
            break;
	case 'W':
            ret = 10;
            break;
	case 'H':
            ret = 11;
            break;
	case 'K':
            ret = 12;
            break;
	case 'R':
            ret = 13;
            break;
        case 'Q':
            ret = 14;
            break;
	case 'N':
            ret = 15;
            break;
	case 'E':
            ret = 16;
            break;
	case 'D':
            ret = 17;
            break;
	case 'S':
            ret = 18;
            break;
        case 'T':
            ret = 19;
    }
    
//    printf("out funhashSimplified\n");
    return ret;
}

unsigned long int insere(cabeca ** hash, char * elem, int pos, int janela, int rep, int gap, long inserir)
{
    unsigned long int i;
    int tam;
    long antHash = 0;
    
//    printf("insere \n");
    
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir * 19; // desloca 5 bits a esquerda
        i       = inserir + funhashSimplified(elem[janela -1]);

    }
    
    antHash = i - funhashSimplified(elem[0])*(pow(19,janela - 1));
    
    int gap_aux 	= 0;
    if (hash[i] == NULL) { //Nao existem nenhum bloco nesta posicao da tabela hash
        hash[i] 	= (cabeca*) calloc(1, sizeof (cabeca)); //Cria uma estrutura no cabeca
        hash[i]->first 	= aloca(elem, pos, (pos + janela - 1), janela); //Seta o primeiro bloco repetitivo nesta posicao da tabela hash
        hash[i]->last 	= hash[i]->first;
    } else {
        lista * aux = hash[i]->last; //Pegamos o ultimo bloco repetitivo da lista
        if ((pos - aux->posicaofim) >= 1 && (pos - aux->posicaofim) <= (1+gap)) { //Fragmento atual e sequencia no bloco repetitivo
	      gap_aux = pos - aux->posicaofim - 1;
	      aux->gaps += gap_aux;
	      aux->lenght++;
	      tam = aux->lenght;
	      aux->posicaofim += (janela+gap_aux); //Aumento o tamanho do bloco        alterado by robson
	    
        } else if ((pos - aux->posicaofim) > (1+gap)) {//Fragmento atual nao sequencia no bloco repetitivo
            if (aux->lenght >= rep) { //Inicia-se um novo bloco repetitivo
                aux->pro = aloca(elem, pos, (pos + janela - 1), janela);
                aux->pro->ant = aux;
                hash[i]->last = aux->pro;
            } else { //Remove o bloco repetitivo atual, pois ele nao tem o tamanho necessirio
                aux->posicaoini = pos;
                aux->posicaofim = pos + janela - 1;
		aux->lenght 	= 1; // Zera o contador do tamanho da repetição
		aux->gaps 	= 0;

            }
            
        }
    }
//    printf("out insere\n");
    return antHash;
}
unsigned long int insereDegen(head ** hashDegen,char * elem , int janela, long inserir)
{
    unsigned long int i;   
    long antHash = 0;
    
//    printf("insereDegen \n");
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir * 19; // desloca 5 bits a esquerda
        i       = inserir + funhashSimplified(elem[janela -1]);

    }
    
    antHash = i - funhashSimplified(elem[0])*(pow(19,janela - 1));
   
//    printf("%lu 506\n",i);
    if(hashDegen[i] == NULL){ // Não existe nunhum bloco nesta posição da tabela hash
      
        hashDegen[i]            = (head*) malloc(sizeof(head));; // cria uma estrutura 
        hashDegen[i]->first      = alocaDegen(elem,i,janela);
        hashDegen[i]->contVariation = 0 ; // Variavel que vai armazenas a qtidade de variações geradas;
        hashDegen[i]->last      = hashDegen[i]->first;
    }
    
    //printf("out insereDegen\n");
    return antHash;
}

void insereDegenVariation(head ** hashDegen, char * elem, int janela, int position)
{
    int i=0;
    unsigned long int hashValue = funhash(elem,janela);
    
 //   printf("insereDegenVariation \n");
    if(hashDegen[hashValue]!= NULL){
        listaDegen * aux = hashDegen[position]->last;
        aux->proximo = alocaDegen(elem,hashValue,janela);
        aux->proximo->anterior = aux;
        hashDegen[position]->contVariation++;
        hashDegen[position]->last = aux->proximo;
                
    } 
    
//    printf("out insereDegenVariation\n");

}


void insereDegeneration(cabeca ** hash, head ** hashDegen, char * elem, int pos, int janela, int rep, int gap, long inserir)
{
    unsigned long int i;
    int cont;
    listaDegen * auxHashDegen;
    lista * auxHash;
    int gap_aux 	= 0;
    
    
 //   printf("insereDegeneration \n");
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir * 19; // desloca 5 bits a esquerda
        i       = inserir + funhashSimplified(elem[janela -1]);

    }
    
    
 //   printf("linha 539 %s %lu \n",elem,i);
    if(hashDegen[i] != NULL){
        auxHashDegen = hashDegen[i]->first;
//	printf("linha 542 \n");
        if(auxHashDegen!=NULL){ // pula primeiro elemento
            auxHashDegen = auxHashDegen->proximo;
        }
        while(auxHashDegen !=NULL){
            if(hash[auxHashDegen->hashvalue]!=NULL){
                auxHash = hash[auxHashDegen->hashvalue]->last;        
                if ((pos - auxHash->posicaofim) >= 1 && (pos - auxHash->posicaofim) <= (1+gap)) { //Fragmento atual e sequencia no bloco repetitivo

                    gap_aux = pos - auxHash->posicaofim - 1;                    
                    auxHash->gaps += gap_aux;
                    auxHash->deg="*";
                    auxHash->posicaofim += (janela+gap_aux); //Aumento o tamanho do bloco        alterado by robson
                    
                }
            }
            
            auxHashDegen = auxHashDegen->proximo;
        }
    }
    
//    printf("out insereDegeneration\n");
}

void generate(head ** hashDegen,char * elem, int iniPos, int janela, int txDegen, int deep, int position){
    char auxSSR[6]; // 6 é o tamanho maximo da janela
    char nucleotideos[] = "GPAVLIMCFYWHKRQNEDST";
    int cont = 0, aux;
    int contProtein = 0;
    
//    printf("generate \n");
    
    
    if(deep < txDegen){
        deep++;
        for(aux = iniPos; aux < janela; aux++){
            strcpy(auxSSR,elem);
	    
	    for(contProtein = 0; contProtein < 20 ; contProtein++){
		auxSSR[aux] = nucleotideos[contProtein];
		if(strcmp(auxSSR,elem)){
		    insereDegenVariation(hashDegen,auxSSR,janela, position);
		    if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}		  
		}
	    }
	    /*
            auxSSR[aux] = nucleotideos[0];
            if(strcmp(auxSSR,elem)){
                 insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }

            auxSSR[aux] = nucleotideos[1];
            if(strcmp(auxSSR,elem)){
                 insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }

            auxSSR[aux] = nucleotideos[2];
            if(strcmp(auxSSR,elem)){
                insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }

            auxSSR[aux] = nucleotideos[3];
            if(strcmp(auxSSR,elem)){
                  insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }
            */
        }
    }
    
//    printf("out generate\n");
    
}
void libera(cabeca ** hash, unsigned long int tam)
{
    int i;
    
//    printf("libera \n");
    for (i = 0; i < tam; i++) { //Percorre toda tabela hash
        if (hash[i] != NULL) {
            lista * aux, *temp;
            aux = hash[i]->first; //Pega o primeiro elemento da lista
            while (aux != NULL) { //Percorre toda lista em cada posicao da tabela hash
                temp = aux;
                aux = aux->pro;
                free(temp);
            }
            free(hash[i]); //Libera memoria do no cabeca
            hash[i] = NULL;
        }
    }
//    printf("out libera\n");
}

void liberaHashDegen(head ** hash, unsigned long int tam)
{
    unsigned long int i;
//    printf("liberaHashDegen \n");
    for (i = 0; i < tam; i++) { //Percorre toda tabela hash
        if (hash[i] != NULL) {
            listaDegen * aux, *temp;
            aux = hash[i]->first; //Pega o primeiro elemento da lista
            while (aux != NULL) { //Percorre toda lista em cada posicao da tabela hash
                temp = aux;
                aux = aux->proximo;
                free(temp);
            }
            free(hash[i]); //Libera memoria do no cabeca
            hash[i] = NULL;
        }
    }
//    printf("out liberaHashDegen\n");
}

int nulo(char * elem, int janela)
{
    int i;
 //   printf("nulo \n");
    for (i = 0; i < janela; i++) { //Percorre todo fragmento
        if (elem[i] == 'X')
            return 0; //Existe pelo menos um 'N' dentro do fragmento
    }
    return 1; //Nao existe 'N' dentro do fragmento
}

__inline int nuloChar(char * elem, int janela){
    
    if(elem[janela - 1]=='X'){
        return 0;
    }
//    printf("out nulo\n");
    return 1;
}

void resultado(FILE * fp, cabeca ** hash, char * nome, int tam, int janela, int rep)
{
    int f, g, i;
    unsigned long int p = pow(20, janela);
    lista * aux, * temp;
    int alread_printed_head = 1;

    f = 1;
    g = 0;
    
 //   printf("resultado \n");
    for (i = 0; i < p; i++) { //Percorre toda tabela hash
        if (hash[i] != NULL) { //Encontra posicao que contem elementos
            aux = hash[i]->first; //Pega o primeiro elemento da lista
            while (aux != NULL) { //Percorre todos bloco repetitivos nesta lista

                int l = aux->lenght;

		if (l >= rep) {
                    if (f) {
                        f = 0;
                    }

                    g = 1;
                    temp = aux; //Guarda primeiro elemento do bloco repetitivo
                    fprintf(fp, "\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s%s",nome, tam, janela, aux->lenght, aux->posicaoini+1, aux->posicaofim+1, aux->gaps, aux->deg,aux->frag);
                }
                aux = aux->pro; //Avanca para o proximo bloco repetitivo
            }
            if (g) {

                g = 0;
            }

        }
    }
    if (f) {
    }
//    printf("out resultado\n");
}
