#include <stdio.h>

#include <omp.h>

int main(int argc, char const *argv[])
{
    #pragma omp parallel                   
    {
        printf("Hello World... from thread = %d\n", omp_get_thread_num());

        #pragma omp barrier

        #pragma omp single
        {
            printf("There are in total %d threads\n", omp_get_num_threads());
        }
    }  

    return 0;
}
