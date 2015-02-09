/**
new Dummy Scheduler with multiple managers
**/
/*
3 functions:
master - manages the managers and sends them them a set of jobs. (size of the set should be around 1000)
manager - Receives the set of jobs and then send them one by one to the workers.
worker - Receives a job from the manager and perform that job.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>

#define POOLSIZE 1024
#define PIVOT 512

MPI_Comm new_comm;

int main(int argc , char * argv[])
{
	int job_start = 0;
	int job_end = 100000;
    /* deifining variables */

    int mpitask_id,number_of_processors;

	/* Intialize MPI*/

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpitask_id);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);

    int color;

    int noProc = number_of_processors - (number_of_processors%POOLSIZE);

    if(mpitask_id == 0)
    {
        color = 0;
    }
    else if ( mpitask_id > noProc)//Assigning colors according to the number of processors available
    {
        if((number_of_processors-noProc)<PIVOT)
        {
            color = ( (mpitask_id - noProc)% (number_of_processors/POOLSIZE)) + 1;
        }
        else color = (number_of_processors/POOLSIZE) + 1;
    }
    else color = (mpitask_id/POOLSIZE)+1;

	//printf("%d       %d\n",mpitask_id, color);

    int new_rank, new_size;

    MPI_Comm_split(MPI_COMM_WORLD,color,mpitask_id, &new_comm); //Creating groups on the basis of color. One of the groups would contain only one processor and this would become the master
    MPI_Comm_rank(new_comm, &new_rank);

    if(mpitask_id == 0)
    {
		int number_of_managers = number_of_processors/POOLSIZE;;

/*Calling the master function with the number of managers depending on the if (total number of processors in MPI_COMM_WORLD) %(mod) POLLSIZE is greater than PIVOT*/
		if(noProc<PIVOT)
		{
			master(number_of_managers,job_start,job_end);
		}
		else
		{
			++number_of_managers;
			master(number_of_managers,job_start,job_end);
		}
    }
    else if(new_rank == 0)
    {
        manager();
    }
    else
    {
        worker();
    }
    MPI_Finalize();
    return 0;
}

void master(int number_of_managers, int job_start, int job_end)
{
    int i,job_sent,job_start=0,job_end=10000; //Set the job size here
    job_sent = 0;
    int job[2],comm_size;
    MPI_Status status;
    /*Send a set of jobs to be performed to all the managers*/
    for(i=0;i<number_of_managers;++i)
    {
        job[0] = job_sent;
        MPI_Recv(&comm_size,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
        job_sent += (comm_size-1);
        job[1] = job_sent;
        int dest = status.MPI_SOURCE;
        MPI_Send(job,2,MPI_INT,dest,1,MPI_COMM_WORLD);
    }
    while(job_sent<job_end)
    {
        MPI_Recv(&comm_size,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);//Wait for a manager to request for jobs
        int dest = status.MPI_SOURCE;
        job[0] = job_sent;
        job_sent += 1+((job_end-job_sent)/comm_size);///I think this would be good since the number of jobs sent would decrease according to the number of jobs left in the queue.
        job[1] = job_sent;
        MPI_Send(job,2,MPI_INT,dest,1,MPI_COMM_WORLD);//Sends job to the manager
    }
    for(i=0;i<number_of_managers;++i)
    {
        MPI_Recv(&comm_size,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);//Wait for a manager to request for jobs
        int dest = status.MPI_SOURCE;
        MPI_Send(job,2,MPI_INT,dest,0,MPI_COMM_WORLD);//Send message to the manager that no more jobs are left.
    }
}

void manager()
{
    MPI_Status status;
    int comm_size,temp,comm_rank,i;
    MPI_Comm_rank(new_comm, &comm_rank);
    MPI_Comm_size(new_comm,&comm_size);

    int job_start,job_end,tag,job_curr;
    int job[2];
    MPI_Send(&comm_size,1,MPI_INT,0,0,MPI_COMM_WORLD);//Requests the master to send jobs
    MPI_Recv(job,2,MPI_INT,0,1,MPI_COMM_WORLD,&status);//Receives a set of jobs to perform from the master
    job_start = job[0];
    job_end = job[1]-1;
    job_curr = job_start;
    for (i=1; i<comm_size; i++)
    {
        MPI_Send(&job_curr, 3, MPI_INT, i, 1, new_comm);//Send a job to the worker
        job_curr++;
    }
    tag = status.MPI_TAG;
    int flag = 0; //Using this flag to ensure that the manager requests master for the job only when a worker is free.
    int dest;
    while(tag)//This loop runs until the master doesn't send the signal to stop.
    {
        while(job_curr<job_end)//Sends the jobs in its queue to the worker who requests for a new job until the queue is not empty.
        {
            if(flag == 1)
            {
                MPI_Recv(&temp,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,new_comm,&status);
                dest = status.MPI_SOURCE;
            }
            flag = 1;
            MPI_Send(&job_curr, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            job_curr++;
        }
        flag = 0;
        MPI_Recv(&temp,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,new_comm,&status);//Receives requests from a worker for work
        dest = status.MPI_SOURCE;
        MPI_Send(&comm_size,1,MPI_INT,0,1,MPI_COMM_WORLD);//Requests master for more jobs
        MPI_Recv(job,2,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);//Receives a set of jobs from the master.
        tag = status.MPI_TAG;
        job_start = job[0];
        job_end = job[1]-1;
        job_curr = job_start;
    }
    /*Send message to the workers signalling 'No More Jobs Left'*/
    MPI_Send(&job_curr, 1, MPI_INT, dest,0, MPI_COMM_WORLD);
    for(i=1;i<(comm_size-1);++i)
    {
        MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, new_comm, &status);
        who = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        MPI_Send(&job, 1, MPI_INT, who, 0, new_comm);
    }
}

void writer()
{
    MPI_Status status;
    int comm_size,temp,comm_rank,i;
    int job;
    MPI_Recv(&job,1,MPI_INT,0,MPI_ANY_TAG,new_comm,&status); //Receives jobs from the manager
    int tag = status.MPI_TAG;
    while(tag)//This loop runs until the manager doesn't send the signal to stop.
    {
        printf("%d",job);///Do Work Here

        MPI_Send(&job,1,MPI_INT,0,1,new_comm);//Request manager for more work
        MPI_Recv(&job,1,MPI_INT,0,MPI_ANY_TAG,new_comm,&status);//Receive work from manager
        tag = status.MPI_TAG;
    }
}
