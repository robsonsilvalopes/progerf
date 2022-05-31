/******************************************************************************

data    : 26 de setembro de 2013
Horário : 10:51
Autor   : Robson

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
    int         hashvalue;      // armazena o valor hash
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
int             funhash         (char * frag, int janela);
__inline int    funhashSimplified(char a);
int             insere          (cabeca ** hash, char * elem, int pos, int janela, int rep, int gap, int inserir);
int             insereDegen     (head ** hashDegen,char * elem, int janela, int inserir);
void            insereDegenVariation(head ** hashDegen,char * elem, int janela, int position);
void            insereDegeneration(cabeca ** hash, head ** hashDegen, char * elem, int pos, int janela, int rep, int gap, int inserir);
void            generate        (head ** hashDegen,char * elem, int iniPos, int janela, int txDegen, int deep, int position);
void            libera          (cabeca ** hash, long tam);
int             nulo            (char * elem, int janela);
__inline int    nuloChar(char * elem, int janela);
void            resultado       (FILE * fp, cabeca ** hash, char * nome, int tam, int janela, int rep);
void            liberaHashDegen (head ** hash, long tam);
fpos_t          pos; //armazena a posição do stream que varre o arquivo fpin
char            nucleotideos[] = "ACTG";

int main(int argc, char ** argv)
{
    FILE * fpin, *fpout;
    int count, f, i, janela, rep, gap,flag, flagEntryDegen,antHash,qtNuclChange;
    int countAux;
    long p;
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
            fprintf(stderr, "Paramento inválido.\n");
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
    
    
    p           = pow(4, janela); //Calcula o tamanho da tabela hash
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
 //                       printf ("-------%i \t %s \t %i \n",countAux, hashDegen[countAux]->first->frag,hashDegen[countAux]->first->hashvalue);                        
                        //generateVariation(hashDegen,countAux,janela,1);
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
    lista * novo = (lista*) calloc(1, sizeof (lista)); //Cria uma nova estruta
    novo->frag = (char*) calloc(janela, sizeof (char)); //Cria um vetor dentro da estrutura
    strcpy(novo->frag, elem); //Copia fragmento para dentro da estrura
    novo->posicaoini = posini; //Seta posicao inicial do fragmento
    novo->posicaofim = posfim;
    novo->gaps = 0;
    novo->lenght = 1;
    novo->deg = (char*) calloc(2, sizeof (char)); //Cria um vetor dentro da estrutura
    novo->deg = "";
    return novo;
}


listaDegen * alocaDegen(char * elem, int valorHash, int janela)
{
    listaDegen * novo = (listaDegen*) calloc(1, sizeof (listaDegen)); //cria uma nova estrutura
    novo->frag     = (char*) calloc(janela, sizeof (char)); //Cria um vetor dentro da estrutura
    strcpy(novo->frag, elem);
    novo->hashvalue = valorHash;
    return novo;
    
}
int funhash(char * frag, int janela)
{
    int i, t, ret = 0;
    t = janela - 1;
    for (i = 0; i < janela; i++) {
        switch (frag[i]) {
            case 'C':
                ret += (pow(2, ((t - i) * 2)));
                break;
            case 'G':
                ret += (pow(2, (((t - i) * 2) + 1)));
                break;
            case 'T':
                ret += (pow(2, (((t - i) * 2) + 1)));
                ret += (pow(2, ((t - i) * 2)));
        }
    }
    return ret;
}

__inline int funhashSimplified(char a){
    int ret = 0;
    switch(a){
        case 'C':
            ret = 1;
            break;
        case 'G':
            ret = 2;
            break;
        case 'T':
            ret = 3;
    }
    return ret;
}

int insere(cabeca ** hash, char * elem, int pos, int janela, int rep, int gap, int inserir)
{
    int i , tam;
    
    
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir << 2;
        inserir = LAST(inserir,janela*2);
        i       = inserir + funhashSimplified(elem[janela -1]);
  //      i = funhash(elem, janela);
    }
    
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
 //           if (((aux->posicaofim - aux->posicaoini + 1) / janela) >= rep) { //Inicia-se um novo bloco repetitivo
	   if ( aux->lenght >= rep) { //Inicia-se um novo bloco repetitivo
                aux->pro = aloca(elem, pos, (pos + janela - 1), janela);
                aux->pro->ant = aux;
                hash[i]->last = aux->pro;
            } else { //Remove o bloco repetitivo atual, pois ele nao tem o tamanho necessirio
                aux->posicaoini = pos;
                aux->posicaofim = pos + janela - 1;
		aux->lenght = 1; // Zera o contador do tamanho da repetição
		aux->gaps = 0;

            }
            
        }
    }
    return i;
}
int insereDegen(head ** hashDegen,char * elem , int janela, int inserir)
{
    int i;    
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir << 2;
        inserir = LAST(inserir,janela*2);
        i       = inserir + funhashSimplified(elem[janela -1]);
   }
    
    if(hashDegen[i] == NULL){ // Não existe nunhum bloco nesta posição da tabela hash
        hashDegen[i]            = (head*) malloc(sizeof(head));; // cria uma estrutura 
        hashDegen[i]->first      = alocaDegen(elem,i,janela);
        hashDegen[i]->contVariation = 0 ; // Variavel que vai armazenas a qtidade de variações geradas;
        hashDegen[i]->last      = hashDegen[i]->first;
    }
    
    return i;
}

void insereDegenVariation(head ** hashDegen, char * elem, int janela, int position)
{
    int i=0;
    int hashValue = funhash(elem,janela);
    if(hashDegen[hashValue]!= NULL){
        listaDegen * aux = hashDegen[position]->last;
        aux->proximo = alocaDegen(elem,hashValue,janela);
        aux->proximo->anterior = aux;
        hashDegen[position]->contVariation++;
        hashDegen[position]->last = aux->proximo;
                
    } 

}


void insereDegeneration(cabeca ** hash, head ** hashDegen, char * elem, int pos, int janela, int rep, int gap, int inserir)
{
    int i,cont;
    listaDegen * auxHashDegen;
    lista * auxHash;
    int gap_aux 	= 0;
    if(inserir < 0){
        i = funhash(elem, janela);
    }else{
        inserir = inserir << 2;
        inserir = LAST(inserir,janela*2);
        i       = inserir + funhashSimplified(elem[janela -1]);
  //      i = funhash(elem, janela);
    }
    if(hashDegen[i] != NULL){
        auxHashDegen = hashDegen[i]->first;
        if(auxHashDegen!=NULL){ // pula primeiro elemento
            auxHashDegen = auxHashDegen->proximo;
        }
        while(auxHashDegen !=NULL){
            if(hash[auxHashDegen->hashvalue]!=NULL){
                auxHash = hash[auxHashDegen->hashvalue]->last;        
                if ((pos - auxHash->posicaofim) >= 1 && (pos - auxHash->posicaofim) <= (1+gap)) { //Fragmento atual e sequencia no bloco repetitivo

                    gap_aux = pos - auxHash->posicaofim - 1;                    
                    auxHash->gaps += gap_aux;
		    auxHash->lenght++;
                    auxHash->deg="*";
                    auxHash->posicaofim += (janela+gap_aux); //Aumento o tamanho do bloco        alterado by robson
                    
                }
            }
            
            auxHashDegen = auxHashDegen->proximo;
        }
    }
}

void generate(head ** hashDegen,char * elem, int iniPos, int janela, int txDegen, int deep, int position){
    char auxSSR[6];
    char  nucleotideos[] = "ACGT";
    int cont = 0, aux;
    
    
    if(deep < txDegen){
        deep++;
        for(aux = iniPos; aux < janela; aux++){
            strcpy(auxSSR,elem);
            auxSSR[aux] = nucleotideos[0];
            if(strcmp(auxSSR,elem)){
 //               printf("%i - %s\n",cont++,auxSSR);
                insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela , txDegen, deep, position);}
            }

            auxSSR[aux] = nucleotideos[1];
            if(strcmp(auxSSR,elem)){
 //               printf("%i - %s\n",cont++,auxSSR);
                insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }

            auxSSR[aux] = nucleotideos[2];
            if(strcmp(auxSSR,elem)){
//                printf("%i - %s\n",cont++,auxSSR);
                insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }


            auxSSR[aux] = nucleotideos[3];
            if(strcmp(auxSSR,elem)){
  //              printf("%i - %s\n",cont++,auxSSR);
                insereDegenVariation(hashDegen,auxSSR,janela, position);
                if(deep < txDegen){ generate(hashDegen,auxSSR, aux+1, janela, txDegen, deep, position);}
            }
        }
    }
    
}
void libera(cabeca ** hash, long tam)
{
    int i;
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
}

void liberaHashDegen(head ** hash, long tam)
{
    int i;
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
}

int nulo(char * elem, int janela)
{
    int i;
    for (i = 0; i < janela; i++) { //Percorre todo fragmento
        if (elem[i] == 'N')
            return 0; //Existe pelo menos um 'N' dentro do fragmento
    }
    return 1; //Nao existe 'N' dentro do fragmento
}

__inline int nuloChar(char * elem, int janela){
    
    if(elem[janela - 1]=='N'){
        return 0;
    }
    return 1;
}

void resultado(FILE * fp, cabeca ** hash, char * nome, int tam, int janela, int rep)
{
    int f, g, i;
    long p = pow(4, janela);
    lista * aux, * temp;
    int alread_printed_head = 1;
//    fprintf(fp, "Id: %s. Sequence Length: %d. Subsequence Length: %d.\n", nome, tam, janela, rep);
    f = 1;
    g = 0;
    for (i = 0; i < p; i++) { //Percorre toda tabela hash
        if (hash[i] != NULL) { //Encontra posicao que contem elementos
            aux = hash[i]->first; //Pega o primeiro elemento da lista
            while (aux != NULL) { //Percorre todos bloco repetitivos nesta lista
//                int k = (aux->posicaofim - aux->posicaoini + 1) % janela;
                int l = aux->lenght;
//                if (!k && (l >= rep)) {
		if (l >= rep) {
                    if (f) {
//			if(alread_printed_head){
//			  fprintf(fp, "Id| Sequence Length| Subsequence Length| Repeats| Initial Block Pos| Final block Pos| Gaps | Subsequence:"); //Imprime cabecalho
//			  alread_printed_head = 0;
//			}
                        f = 0;
                    }

                    g = 1;
                    temp = aux; //Guarda primeiro elemento do bloco repetitivo
                    fprintf(fp, "\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s%s",nome, tam, janela, aux->lenght, aux->posicaoini+1, aux->posicaofim+1, aux->gaps, aux->deg,aux->frag);
/*   
		    aux->posicaoini = aux->posicaofim - aux->posicaoini + 1;
                    while (aux->posicaoini > 0) {
                        fprintf(fp, "%s", temp->frag); //Imprime o ultimo fragmento
                        aux->posicaoini = aux->posicaoini - janela;
                    }
*/
                }
                aux = aux->pro; //Avanca para o proximo bloco repetitivo
            }
            if (g) {
//                fprintf(fp, "\n");
                g = 0;
            }

        }
    }
    if (f) {
 //       fprintf(fp, "Fragment size: %d.\n", janela);
 //       fprintf(fp, "There are not results for these search parameters.\n");
    }
}
