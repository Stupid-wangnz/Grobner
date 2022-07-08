#include "Grobner_Matrix.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<pthread.h>
#include<windows.h>
#include<sys/time.h>
#include "omp.h"
#include <mpi.h>

using namespace std;

typedef struct {
    int t_id;
    int n;
    vector<int> v;
}threadParam_t;

int NUM_THREADS=4;

int n=8399;
int eliminator_count=106,elminated_element=4535;

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
int read_eliminator(fstream&infile_eliminator){
    if(!infile_eliminator.is_open()){
        cout<<"err open";
        return -1;
    }
    /*
     * 每次读入五行消元子进行消元
     *
     * @return 返回消元子的数量
     */

    eliminator.clear();
    string line;
    int num=0;
    while(num<5&&getline(infile_eliminator,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminator.input_line(num,input_line);
        num++;
    }

    return num;
}
void func(){
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
        return ;
    }
    vector<int>eliminator_row_index(n,-1);
    vector<int>eliminated_element_row_index(elminated_element,-1);

    //read the record of eliminated element
    string line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminatedElement.input_line(record,input_line);
        eliminated_element_row_index[record]=record;
        record++;
    }

    cout<<"the eliminated element is:"<<endl;
    for(int i=0;i<eliminated_element_row_index.size();i++){
        cout<<eliminated_element_row_index[i]<<" ";
    }
    cout<<endl;
    //read the record of eliminator

    record=0;
    line="";
    while(getline(infile_eliminator,line)){
        stringstream ss(line);
        int temp;
        ss>>temp;
        eliminator_row_index[temp]=temp;
    }
    cout<<"the eliminator is:"<<endl;
    for(int i=0;i<eliminator_row_index.size();i++){
        cout<<eliminator_row_index[i]<<" ";
    }
    cout<<endl;

    infile_eliminated_element.close();
    infile_eliminator.close();
    //串行的版本
    /*long long head,tail,freq;
    double time=0;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);*/

    infile_eliminator.open(R"(F:\CL_WorkSpace\Grobner\input1.txt)");

    cout<<"begin to eliminate"<<endl;
    //消元计算
    while(!eliminated_element_row_index.empty()){
        int num=read_eliminator(infile_eliminator);
        if(num==0||num==-1){
            break;
        }
        //消元子消元
        for(int i=0;i<num;i++){
            int row_index=eliminator.row_index[i];
            int eliminator_=eliminator.get_max_bit(row_index);
            for(int j=0;j<eliminated_element_row_index.size();j++){
                if(eliminated_element_row_index[j]==-1){
                    continue;
                }
                if(eliminator_==eliminatedElement.get_max_bit(eliminated_element_row_index[j])){
                    //消元子消元
                    int new_max_bit=eliminatedElement.xor_line(eliminator,i,j);
                    if(new_max_bit==0){
                        //被完全消元
                        eliminated_element_row_index[j]=-1;
                        continue;
                    }
                    if(eliminator_row_index[new_max_bit]==-1){
                        //升格成为新的消元子
                        cout<<"new eliminator:"<<new_max_bit<<endl;
                        eliminatedElement.print_line(j);
                        eliminator_row_index[new_max_bit]=new_max_bit;
                        eliminated_element_row_index[j]=-1;
                    }
                }
            }
        }
        //更新需要消元的行
        vector<int>new_eliminated_element_row_index;
        for(int i=0;i<eliminated_element_row_index.size();i++){
            if(eliminated_element_row_index[i]!=-1){
                new_eliminated_element_row_index.push_back(eliminated_element_row_index[i]);
            }
        }
        eliminated_element_row_index=new_eliminated_element_row_index;
        if(eliminated_element_row_index.empty()){
            break;
        }
    }
    infile_eliminator.close();
    /*QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time+=(tail-head)*1000.0/freq;
    cout<<"time:"<<time<<endl;*/

    //unc_pthread();

    return ;
}
int Serial(){
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

    //func_pthread();

    return 0;
}
int OpenMP(){
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


    long long head,tail,freq;
    double time=0;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    while(eliminatedElement.row_index.size()>0)
    {
        //仍有被消元子没被消元完
        vector<int> row_index=eliminatedElement.row_index;
        int i,old_max_bit,new_max_bit;
        #pragma omp parallel for private(i,old_max_bit,new_max_bit)
        for(i=0;i<row_index.size();i+=1){
            old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
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
        //线程计算完毕，更新eliminator
        for(i=0;i<row_index.size();i++){
            new_max_bit=eliminatedElement.get_max_bit(row_index[i]);
            if(new_max_bit==0){
                //消元子消元完毕
                eliminatedElement.row_index[i]=-1;
                continue;
            }
            if(eliminator.row_index[new_max_bit]==-1){
                eliminator.row_index[new_max_bit]=new_max_bit;
                cout<<"new eliminator:"<<endl;
                eliminatedElement.print_line(row_index[i]);
                eliminatedElement.row_index[i]=-1;
                for(int k=0;k<eliminatedElement.m_;k++){
                    eliminator.matrix[new_max_bit][k]=eliminatedElement.matrix[row_index[i]][k];
                }
            }

        }
        vector<int>new_row_index;
        for(i=0;i<eliminatedElement.row_index.size();i++){
            if(eliminatedElement.row_index[i]!=-1){
                new_row_index.push_back(eliminatedElement.row_index[i]);
            }
        }
        eliminatedElement.row_index=new_row_index;
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<"time:"<<time<<endl;


    return 0;
}
void readdata(){
    ifstream infile_eliminator;
    infile_eliminator.open("F:\\CL_WorkSpace\\test\\7_8399_6375_4535\\1.txt");
    ifstream infile_eliminated_element;
    infile_eliminated_element.open("F:\\CL_WorkSpace\\test\\7_8399_6375_4535\\2.txt");

    if(!infile_eliminator.is_open()||!infile_eliminated_element.is_open()){
        cout<<"err open";
        return ;
    }
    string line="";
    while(getline(infile_eliminator,line)){
        istringstream ss(line);
        int temp;
        int index=-1;
        while(ss>>temp){
            if(index==-1)
                index=temp;
            eliminator.set_bit(index,temp);
        }
        eliminator.row_index[index]=index;
    }
    infile_eliminator.close();
    line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        istringstream ss(line);
        int temp;
        while(ss>>temp){
            eliminatedElement.set_bit(record,temp);
        }
        eliminatedElement.row_index[record]=record;
        record++;
    }
    infile_eliminated_element.close();
}
void super(int rank,int num_proc)
{
    vector<int> row_index=eliminatedElement.row_index;
    int i,old_max_bit,new_max_bit=0;
//#pragma omp parallel for private(i,old_max_bit,new_max_bit)
    for(i=rank;i<row_index.size();i+=num_proc){
        old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
        while(old_max_bit>=0&&eliminator.row_index[old_max_bit]!=-1){
            //消元子还没消元完
            new_max_bit=eliminatedElement.xor_line(eliminator,old_max_bit,row_index[i]);
            old_max_bit=new_max_bit;
            //如果new_max_bit不在eliminator的row_index中，则被消元子消元完毕
            if(new_max_bit<0){
                break;
            }

        }
    }

}
vector<int>new_eliner;
void MPIfunc(){
    int rank;
    int num_proc;//进程数
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//获取当前进程号
    cout<<rank<<endl;

    if(rank==0){
        //0进程分发数据
        timeval t_start;
        timeval t_end;
        gettimeofday(&t_start, NULL);
        int sign;
        do{
            /*cout<<"---------before----------"<<endl;
            for(int i=0;i<eliminatedElement.row_index.size();i++){
                eliminatedElement.print_line(eliminatedElement.row_index[i]);
            }
            cout<<"--------------------------------"<<endl;*/

            for(int i=0;i<eliminatedElement.row_index.size();i++){
                int flag=i%num_proc;
                if(flag==rank){
                    continue;
                }
                else
                    MPI_Send(eliminatedElement.matrix[eliminatedElement.row_index[i]],eliminatedElement.m_,MPI_INT,flag,0,MPI_COMM_WORLD);
            }
            //消元
            super(rank,num_proc);
            cout<<0<<endl;
            //接收新的被消元子
            for (int i = 0; i <eliminatedElement.row_index.size(); i++)
            {
                int flag = i % num_proc;

                if (flag == rank)
                    continue;
                else

                    MPI_Recv(eliminatedElement.matrix[eliminatedElement.row_index[i]], eliminatedElement.m_, MPI_INT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            /*cout<<"---------after----------"<<endl;
            for(int i=0;i<eliminatedElement.row_index.size();i++){
                eliminatedElement.print_line(eliminatedElement.row_index[i]);
            }
            cout<<"--------------------------------"<<endl;*/
            //进程计算完毕，更新eliminator
            vector<int> row_index=eliminatedElement.row_index;
            int i,new_max_bit;
            for(i=0;i<row_index.size();i++){
                new_max_bit=eliminatedElement.get_max_bit(row_index[i]);
                if(new_max_bit<0){
                    //消元子消元完
                    /*cout<<"kong"<<endl;*/
                    eliminatedElement.row_index[i]=-1;
                    continue;
                }
                else {
                    if (eliminator.row_index[new_max_bit] == -1) {
                        eliminator.row_index[new_max_bit] = new_max_bit;
                        new_eliner.push_back(new_max_bit);
                        /*cout << "new eliminator: "<<row_index[i] << endl;
                        eliminatedElement.print_line(row_index[i]);*/
                        eliminatedElement.row_index[i] = -1;
                        for (int k = 0; k < eliminatedElement.m_; k++) {
                            eliminator.matrix[new_max_bit][k] = eliminatedElement.matrix[row_index[i]][k];
                        }
                    }
                    else
                        break;
                }
            }
            vector<int>new_row_index;
            for(i=0;i<eliminatedElement.row_index.size();i++){
                if(eliminatedElement.row_index[i]!=-1){
                    new_row_index.push_back(eliminatedElement.row_index[i]);
                }
            }
            eliminatedElement.row_index=new_row_index;

            sign=new_row_index.size()==0?0:1;
            //向其他进程更新新的消元子
            for(int f=1;f<num_proc;f++){
                int t=new_row_index.size();
                MPI_Send(&t,1,MPI_INT,f,2,MPI_COMM_WORLD);
                MPI_Send(&new_row_index[0],int(new_row_index.size()),MPI_INT,f,3,MPI_COMM_WORLD);
                int size=new_eliner.size();
                MPI_Send(&size,1,MPI_INT,f,4,MPI_COMM_WORLD);
                MPI_Send(&new_eliner[0],size,MPI_INT,f,5,MPI_COMM_WORLD);
                for(int h=0;h<size;h++){
                    MPI_Send(eliminator.matrix[new_eliner[h]],eliminator.m_,MPI_INT,f,6,MPI_COMM_WORLD);
                }
            }
            new_eliner.clear();
            num_proc=min(num_proc,(int)eliminatedElement.row_index.size());
        } while (sign==1);

        gettimeofday(&t_end, NULL);
        cout << "super time cost: "
             << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;

    }
    else{
        //其余线程先接收分发的数据

        int sign;
        do{
            for (int i = rank; i <eliminatedElement.row_index.size() ; i += num_proc)
            {
                MPI_Recv(eliminatedElement.matrix[eliminatedElement.row_index[i]], eliminatedElement.m_, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            super(rank,num_proc);
            cout<<rank<<endl;
            for (int i = rank; i < eliminatedElement.row_index.size(); i += num_proc)
            {
                MPI_Send(eliminatedElement.matrix[eliminatedElement.row_index[i]], eliminatedElement.m_, MPI_INT, 0, 1, MPI_COMM_WORLD);
            }

            //接收新的消元子和更新的列索引
            int new_row_index_size=0;
            MPI_Recv(&new_row_index_size,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            eliminatedElement.row_index.resize(new_row_index_size);
            MPI_Recv(&eliminatedElement.row_index[0],new_row_index_size,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            int new_eliner_size=0;
            MPI_Recv(&new_eliner_size,1,MPI_INT,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            new_eliner.resize(new_eliner_size);
            MPI_Recv(&new_eliner[0],new_eliner_size,MPI_INT,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for(int t=0;t<new_eliner_size;t++){
                eliminator.row_index[new_eliner[t]]=new_eliner[t];
                MPI_Recv(eliminator.matrix[new_eliner[t]],eliminator.m_,MPI_INT,0,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            new_eliner.clear();
            num_proc=min(num_proc,(int)eliminatedElement.row_index.size());
            if(rank>=num_proc)
                return;
            sign=new_row_index_size==0?0:1;
        }while(sign==1);

    }


}
int main(){
    MPI_Init(0,0);
    readdata();
    MPIfunc();
    MPI_Finalize();

}