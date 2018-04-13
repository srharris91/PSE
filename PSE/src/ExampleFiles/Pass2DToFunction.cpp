#include<iostream>
using namespace std;

//function to display 2D array
void show(int **arr,int x,int y)
{
    int i,j;
    for(i=0;i<x;i++)
    { 
        for(j=0;j<y;j++)
        {
            cout<<arr[i][j]<<" ";
        }
        cout<<endl;
    }
}

//main program
int main()
{  
    int n,m;
    //input number of rows and columns
    cout<<"Enter No. of rows: ";
    cin>>n;
    cout<<"Enter No. of columns: ";
    cin>>m;

    //pointer to 2D array
    int **A=new int*[n];

    //pointer initialization
    for(int i=0;i<n;i++)
    {
        A[i]=new int[m];
    }

    //input array elements
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            cin>>A[i][j];
        }
    }

    //display 2D array
    show(A,n,m);

    return 0;
}
