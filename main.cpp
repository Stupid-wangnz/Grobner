#include "Grobner_Matrix.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<pthread.h>
#include<ctime>
#include<cstdlib>
#include<windows.h>
#include<time.h>

using namespace std;

typedef struct {
    int t_id;
    int n;
    vector<int> v;
}threadParam_t;


int NUM_THREADS=4;

int n=2362;
int eliminator_count=170,elminated_element=453;

Grobner_Matrix eliminator(n,n);
Grobner_Matrix eliminatedElement(elminated_element,n);


void * threadfunc(void*param){
    threadParam_t * p=(threadParam_t*)param;
    int t_id=p->t_id;
    int step=NUM_THREADS;
    vector<int> row_index=p->v;
    for(int i=t_id;i<row_index.size();i+=step){
        int old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
        int new_max_bit;
        while(eliminator.row_index[old_max_bit]!=-1){
            //消元子还没消元完
            new_max_bit=eliminatedElement.xor_line(eliminator,old_max_bit,row_index[i]);
            //如果new_max_bit不在eliminator的row_index中，则被消元子消元完毕
            if(eliminator.row_index[new_max_bit]==-1){
                break;
            }
            old_max_bit=new_max_bit;
        }
    }

}
void func_pthread(){
    fstream infile_eliminator(R"(F:\CL_WorkSpace\Grobner\input1.txt)");
    fstream infile_eliminated_element(R"(F:\CL_WorkSpace\Grobner\input2.txt)");

    if(!infile_eliminator.is_open()||!infile_eliminated_element.is_open()){
        cout<<"err open";
        return ;
    }

    string line;
    while(getline(infile_eliminator,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminator.input_line(input_line[0],input_line);
    }

    line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        stringstream ss(line);
        vector<int>input_line;

        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }

        eliminatedElement.input_line(record,input_line);
        record++;
    }
    infile_eliminated_element.close();
    infile_eliminator.close();

    threadParam_t *param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
    }
    pthread_t *pthreads=new pthread_t[NUM_THREADS];

    long long head,tail,freq;
    double time=0;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    while(eliminatedElement.row_index.size()>0){
        //仍有被消元子没被消元完
        vector<int> row_index=eliminatedElement.row_index;
        for(int j=0;j<NUM_THREADS;j++){
            param[j].v=row_index;
        }
        for(int j=0;j<NUM_THREADS;j++){
            pthread_create(&pthreads[j],NULL,threadfunc,&param[j]);
        }
        for(int j=0;j<NUM_THREADS;j++){
            pthread_join(pthreads[j],NULL);
        }
        //线程计算完毕，更新eliminator
        for(int i=0;i<row_index.size();i++){
            int new_max_bit=eliminatedElement.get_max_bit(row_index[i]);
            if(new_max_bit==0){
                //消元子消元完毕
                eliminatedElement.row_index[i]=-1;
                continue;
            }
            if(eliminator.row_index[new_max_bit]==-1){
                eliminator.row_index[new_max_bit]=new_max_bit;
                eliminatedElement.row_index[i]=-1;
                for(int k=0;k<eliminatedElement.m_;k++){
                    eliminator.matrix[new_max_bit][k]=eliminatedElement.matrix[row_index[i]][k];
                }
            }
        }
        vector<int>new_row_index;
        for(int i=0;i<eliminatedElement.row_index.size();i++){
            if(eliminatedElement.row_index[i]!=-1){
                new_row_index.push_back(eliminatedElement.row_index[i]);
            }
        }
        eliminatedElement.row_index=new_row_index;
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<"time:"<<time<<endl;

    delete []param;
    delete []pthreads;
}

int main(){

    /*
     * n是总行数
     *
     * eliminator是初始时消元子的数量
     *
     * eliminated_element是初始时被消元子的数量
     */

    fstream infile_eliminator(R"(F:\CL_WorkSpace\Grobner\input1.txt)");
    fstream infile_eliminated_element(R"(F:\CL_WorkSpace\Grobner\input2.txt)");

    if(!infile_eliminator.is_open()||!infile_eliminated_element.is_open()){
        cout<<"err open";
        return 0;
    }

    string line;
    while(getline(infile_eliminator,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminator.input_line(input_line[0],input_line);
    }

    line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        stringstream ss(line);
        vector<int>input_line;

        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }

        eliminatedElement.input_line(record,input_line);
        record++;
    }
    infile_eliminated_element.close();
    infile_eliminator.close();

    //串行的版本
    long long head,tail,freq;
    double time=0;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    cout<<"eliminator:"<<endl;
    int num=0;
    while(num*4<elminated_element){
        //仍有被消元子没被消元完
        vector<int> row_index=eliminatedElement.get_line(num,4);
        for(int i=0;i<row_index.size();i++){
            int old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
            int new_max_bit;
            while(eliminator.row_index[old_max_bit]!=-1){
                //消元子还没消元完
                new_max_bit=eliminatedElement.xor_line(eliminator,old_max_bit,row_index[i]);
                if(new_max_bit==0){
                    //消元子消元完了
                    //cout<<row_index[i]<<": all0"<<endl;
                    break;
                }
                //如果new_max_bit不在eliminator的row_index中，则被消元子消元完毕
                if(eliminator.row_index[new_max_bit]==-1){
                    for(int k=0;k<eliminatedElement.m_;k++){
                        eliminator.matrix[new_max_bit][k]=eliminatedElement.matrix[row_index[i]][k];
                    }
                    eliminator.row_index[new_max_bit]=new_max_bit;
                    break;
                }
                old_max_bit=new_max_bit;
            }
        }
        num++;
    }

    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time+=(tail-head)*1000.0/freq;
    cout<<"time:"<<time<<endl;

    func_pthread();

    return 0;
}

