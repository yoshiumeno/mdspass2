#include <stdio.h>
#include <math.h>

#define n 3             //n-1次の多項式で近似（１以上）
#define N 3             //ガウスの消去法における未知数の数（上のnと同じ値にすること）
#define S 5             //データの個数
#define CHECK 0 //ガウスの消去法における三角行列のチェック用

void gauss(double a[N][N+1],double p[N]);

void sai(double x[S],double y[S])
{
        int i,j,k;
        double X,Y;
        double A[n][n+1],xx[N];

        //FILE *output1;
        //FILE *output2;
        //output1=fopen("output1.data","w");
        //output2=fopen("output2.data","w");

/*初期化*/
        for(i=0;i<n;i++) {
                for(j=0;j<n+1;j++) {
                        A[i][j]=0.0;
                }
        }

/*ガウスの消去法で解く行列の作成*/
        for(i=0;i<n;i++) {
                for(j=0;j<n;j++) {
                        for(k=0;k<S;k++) {
                                A[i][j]+=pow(x[k],i+j);
                        }
                }
        }
        for(i=0;i<n;i++) {
                for(k=0;k<S;k++) {
                        A[i][n]+=pow(x[k],i)*y[k];
                }
        }
/*ガウスの消去法の実行（配列xxは解、すなわち多項式の係数を入れるためのもの）*/
        gauss(A,xx);
	printf("%e %e %e\n",xx[0],xx[1],xx[2]);

/*GNUPLOTで表示するために最小２乗法による関数のデータをファイル保存*/
        for(X=x[0]-10.0;X<x[S]+10.0;X+=0.01) {
                Y=0.0;
                for(i=0;i<N;i++) {
                        Y+=xx[i]*pow(X,i);
                }
                //fprintf(output1,"%f %f\n",X,Y);
        }

/*GNUPLOTで表示するために、最小２乗法に使われたデータを保存*/
        for(i=0;i<S;i++) {
	  //fprintf(output2,"%f %f\n",x[i],y[i]);
        }

        //fclose(output1);
        //fclose(output2);

}

void gauss(double a[N][N+1],double xx[N])
{
        int i,j,k,l,pivot;
        double x[N];
        double p,q,m,b[1][N+1];

        for(i=0;i<N;i++) {
                m=0;
                pivot=i;

                for(l=i;l<N;l++) {
                        if(fabs(a[l][i])>m) {   //i列の中で一番値が大きい行を選ぶ
                                m=fabs(a[l][i]);
                                pivot=l;
                        }
                }

                if(pivot!=i) {                          //pivotがiと違えば、行の入れ替え
                        for(j=0;j<N+1;j++) {
                                b[0][j]=a[i][j];        
                                a[i][j]=a[pivot][j];
                                a[pivot][j]=b[0][j];
                        }
                }
        }

        for(k=0;k<N;k++) {
                p=a[k][k];              //対角要素を保存
                a[k][k]=1;              //対角要素は１になることがわかっているから

                for(j=k+1;j<N+1;j++) {
                        a[k][j]/=p;
                }

                for(i=k+1;i<N;i++) {
                        q=a[i][k];

                        for(j=k+1;j<N+1;j++) {
                                a[i][j]-=q*a[k][j];
                        }
                a[i][k]=0;              //０となることがわかっているところ
                }
        }

//解の計算
        for(i=N-1;i>=0;i--) {
                x[i]=a[i][N];
                for(j=N-1;j>i;j--) {
                        x[i]-=a[i][j]*x[j];
                }
        }

//行列が最後どうなったか見たいときに実行
#if CHECK==1
        for(i=0;i<N;i++) {
                for(j=0;j<N+1;j++) {
                        printf("%10.3f",a[i][j]);
                }
                printf("\n");
                
        }
#endif

        printf("解は\n");
        for(i=0;i<N;i++) {
                printf("%f\n",x[i]);
                xx[i]=x[i];
        }

}


/*
int main(void)
{
        double x[S]={1.0,2.0,3.0,4.0,5.0},y[S]={7.987,2.986,1.998,2.224,5.678};

//データの横軸に値する数値を配列xに小さい順に入れ、
//それぞれの横軸の値における縦軸の値を配列yに入れ関数saiに渡す
        sai(x,y);

        return 0;
}
*/
