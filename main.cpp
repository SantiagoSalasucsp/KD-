//
//  main.cpp
//  FINAL_EDA
//
//  Created by Santiago Salas Sotillo on 18/11/25.
//

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
using namespace std;

 
long long word_claud =10;
#define numero_ladrillos_base 16
#define timepo 2



#define PATH "/Users/santiagosalas/Desktop/python/mensaje.txt"
#define numero_de_lectura 10

template <class T>
struct Cvector{
    T* M;
    T* ini;
    T* fin;
    int tam;
    int r;

    Cvector(int t){
        M=new T[t];
        tam=t;
        ini=M;
        fin=M;
        r=t;
    }

    void push_back(T val){
        if(fin==M+tam-1){
            *fin=val;
            T*tmp=M;
            M=new T[tam*2];
            T* p=M;
            T* q=ini;
            while(q<fin){
                *p=*q;
                q++;
                p++;
            }
            *p=val;
            fin=p+1;
            ini=M;
            tam*=2;
            delete []tmp;

        }
        else{
            *fin=val;
            fin++;
        }
    }


    void push_frond(T val){
        if(fin==M+tam-1){
            T*tmp =M;
            M=new T[tam*2];
            T*p=M+1;
            T* q=ini;
            while(q<fin){
                *p=*q;
                p++;
                q++;
            }
            *M=val;
            ini=M;
            fin=p;
            tam=tam*2;
            delete []tmp;


        }
        else{
            T* p=fin;
            while(p>ini){
                *p=*(p-1);
                p--;
            }
            *ini=val;
            fin++;
        }
    }


    void pop_back(){
        if(fin==M+(tam/2) && ((tam/2)!=r)){
            T*W=M;
            M= new T[tam/2];
            T*p=M;
            T*q=ini;
            while(q<fin-1){
                *p=*q;
                p++;
                q++;

            }
            ini=M;
            fin=p;
            delete[] W;
            tam=tam/2;
        }
        else{
            fin--;
        }

    }


    void pop_frond(){

        if(fin==M+(tam/2) && ((tam/2)>=r)){
            T*W=M;
            M= new T[tam/2];
            T*p=M;
            T*q=ini+1;
            while(q<fin){
                *p=*q;
                q++;
                p++;
            }
            fin=p;
            ini=M;
            tam=tam/2;
            delete[] W;

        }
        else{
            T* p=ini;
            while(p<fin-1){
                *p=*(p+1);
                p++;
            }
            fin--;
        }

    }


};





template<class T>
struct Cola {
    
    int tam;
    T* A;
    T* head;
    T* tail;
    int count;

    Cola(int n) {
        tam = n;
        A = new T[tam];
        
        head = A;
        tail = A;
        
        count = 0;
    }

    bool push(const T& val) {
        if (count == tam){
            return false;
        }
        
        *head = val;
        head++;
        
        if (head == A + tam){
            head = A;
        }
        
        count++;
        return true;
    }

    bool pop(T& val) {
        
        if (count == 0)
            return false;
        
        val = *tail;
        tail++;
        
        if (tail == A + tam)
            tail = A;
        
        count--;
        return true;
    }

    ~Cola() {
        delete[] A;
    }
    
    
};












class par {
public:
    string palabra;
    int frecuencia;

    par(string p, int f) : palabra(p), frecuencia(f) {}

    string enviar_a_python() const {
        return palabra + "," + to_string(frecuencia) + "\n";
    }
    
    
    par& operator=(const par& x) {
        if (this != &x) {
            palabra = x.palabra;
            frecuencia = x.frecuencia;
        }
        return *this;
    }

    
    
};









































































int longitud(char* p) {
    int tam;
    tam=0;
    
    while (p[tam] != '\0'){
        tam++;
    }
    
    return tam;
}



bool sonIguales(const char* str1, const char* str2) {
    int i;
    
    i=0;
    
    while (str1[i] != '\0' && str2[i] != '\0') {
        
        if (str1[i] != str2[i]){
            
            return 0;
        }
        
        i++;
    }
    
    return str1[i] == str2[i];
    
}

void copiar(char* dest, const char* src) {
    int i;
    i=0;
    while (src[i] != '\0') {
        dest[i] = src[i];
        i++;
    }
    dest[i] = '\0';
}

char aMinuscula(char c) {
    
    if (c >= 'A' && c <= 'Z'){
        return c + 32;
    }
    return c;
    
}



bool esLetra(char c) {
    
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

bool esNumero(char c) {
    
    return c >= '0' && c <= '9';
}



bool esStopword(const char* palabra) {
    const char* stopwords[] = {
        "the", "be", "to", "of", "and", "a", "in", "that", "have", "i",
        "it", "for", "not", "on", "with", "he", "as", "you", "do", "at",
        "this", "but", "his", "by", "from", "they", "we", "say", "her", "she",
        "or", "an", "will", "my", "one", "all", "would", "there", "their",
        "what", "so", "up", "out", "if", "about", "who", "get", "which", "go",
        "me", "when", "make", "can", "like", "time", "no", "just", "him", "know",
        "take", "people", "into", "year", "your", "good", "some", "could", "them",
        "see", "other", "than", "then", "now", "look", "only", "come", "its", "over",
        "think", "also", "back", "after", "use", "two", "how", "our", "work",
        "first", "well", "way", "even", "new", "want", "because", "any", "these",
        "give", "day", "most", "us", "is", "was", "are", "been", "has", "had",
        "were", "said", "did", "having", "may", "should", "am", "being", "more"
    };
    
    int numStopwords = 98;
    for (int i = 0; i < numStopwords; i++) {
        if (sonIguales(palabra, stopwords[i])) return true;
    }
    return false;
}


void lematizar(char* palabra) {
    int len = longitud(palabra);
    
    // Regla: -ing -> quitar
    if (len > 4 && palabra[len-3] == 'i' && palabra[len-2] == 'n' && palabra[len-1] == 'g') {
        palabra[len-3] = '\0';
        len -= 3;
        // Si termina en consonante duplicada, quitar una
        if (len >= 2 && palabra[len-1] == palabra[len-2] && !esLetra(palabra[len-1])) {
            palabra[len-1] = '\0';
        }
        return;
    }
    
    // Regla: -ed -> quitar
    if (len > 3 && palabra[len-2] == 'e' && palabra[len-1] == 'd') {
        palabra[len-2] = '\0';
        return;
    }
    
    // Regla: -s o -es (plurales)
    if (len > 3) {
        if (palabra[len-1] == 's') {
            if (len > 4 && palabra[len-2] == 'e') {
                palabra[len-2] = '\0'; // -es
            } else {
                palabra[len-1] = '\0'; // -s
            }
            return;
        }
    }
    
    // Regla: -ly -> quitar
    if (len > 4 && palabra[len-2] == 'l' && palabra[len-1] == 'y') {
        palabra[len-2] = '\0';
        return;
    }
    
    // Regla: -tion -> t
    if (len > 5 && palabra[len-4] == 't' && palabra[len-3] == 'i' &&
        palabra[len-2] == 'o' && palabra[len-1] == 'n') {
        palabra[len-4] = 't';
        palabra[len-3] = '\0';
        return;
    }
}























int main() {
    par datos[] = {
        par("uno", 1),
        par("dos", 2),
        par("tres", 3),
        par("cuatro", 4),
        par("cinco", 5),
        par("seis", 6),
        par("siete", 7),
        par("ocho", 8),
        par("nueve", 9),
        par("diez", 10)
    };

    par datos2[] = {
        par("A", 1),
        par("B", 2),
        par("C", 3),
        par("D", 4),
        par("E", 5),
        par("F", 6),
        par("G", 7),
        par("H", 8),
        par("I", 9),
        par("J", 10)
    };

    ofstream initial_f(PATH, ios::trunc);
    initial_f.close();

    ofstream f(PATH, ios::app);

    for (int i = 0; i < numero_de_lectura; ++i) {
        f << datos[i].enviar_a_python();
        f.flush();
        sleep(timepo);
    }

    string comando;
    cin >> comando;

    if (comando == "clud") {
        f << "###\n";
        f.flush();

        int nuevo_valor;
        cin >> nuevo_valor;
        word_claud = nuevo_valor;

        f << "$$$\n";
        f.flush();

        for (int i = 0; i < numero_de_lectura; ++i) {
            f << datos2[i].enviar_a_python();
            f.flush();
            sleep(timepo);
        }
    }

    while (true) {
        sleep(5);
    }

    return 0;
}

