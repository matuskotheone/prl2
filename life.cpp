/***
 * Author: Matus Gazdik (xgazdi04)
 * Date: 2024-04-23
 * Project: Parralel representation of the game of life using OpenMPI
 * Description: The Program gets 2 Arguments, the first one is the file with the initial state of the game of life, the second one is the number of steps the game should be played.
 * The program prints the end state of the game on the standard output.
 * The size of word is given by the size of the matrix in input file.
 * All the walls are considered as dead cells.
 * Should work with any number of processes. Tested with 4
 ***/


#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <queue>
#include <unistd.h>

#include <mpi.h>

// Structure to represent a row has a vector of integers and the row number
struct row {
    std::vector<int> data;
    int row_number;
};
typedef struct row row_t;

// Function to read the numbers from the file
// Returns a vector of vectors of integers
std::vector<std::vector<int>> read_numbers(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<int>> numbers;
    while(std::getline(file, line)) {
        std::vector<int> number;
        for(int i = 0; i < line.size(); i++) {
            number.push_back(line[i] == '1');
        }
        numbers.push_back(number);
    }
    return numbers;
}

// Prints the row to the standard output
void print_row(row_t row) {
    for(int j = 0; j < row.data.size(); j++) {
        std::cout << (row.data[j] ? "1" : "0");
    }
    std::cout << std::endl;
}

// Prints the rows worked on by current process
void print_rows(std::vector<row_t> rows, int rank) {
    for(int i = 0; i < rows.size(); i++) {
        std::cout << rank << ": " ;
        print_row(rows[i]);
    }
}


// Function gets 3 rows the current row to calculet above row and below row
// Returns the current row after 1 cycle of the game of life
row_t calc_row(row_t curr, row_t up, row_t down) {
    row_t new_row;
    new_row.data.resize(curr.data.size());
    new_row.row_number = curr.row_number;
    for (int i = 0; i < curr.data.size(); i++) {// cycles over the row elements and calculates the new value
        // the count is used to count the number of alive cells around the current cell
        int count = 0;
        if (i == 0) {// if the current cell is on the left wall
            count += up.data[i];
            count += up.data[i + 1];
            count += curr.data[i + 1];
            count += down.data[i];
            count += down.data[i + 1];
        }
        else if (i == curr.data.size() - 1) { // if the current cell is on the right wall
            count += up.data[i - 1];
            count += up.data[i];
            count += curr.data[i - 1];
            count += down.data[i - 1];
            count += down.data[i];
        }
        else {// if it is somewhere in the middle
            count += up.data[i - 1];
            count += up.data[i];
            count += up.data[i + 1];
            count += curr.data[i - 1];
            count += curr.data[i + 1];
            count += down.data[i - 1];
            count += down.data[i];
            count += down.data[i + 1];
        }
        if (curr.data[i] == 1) { // if the current cell is alive
            if (count < 2) {
                new_row.data[i] = 0;
            }
            else if (count == 2 || count == 3) {
                new_row.data[i] = 1;
            }
            else {
                new_row.data[i] = 0;
            }
        }
        else {// if the current cell is dead
            if (count == 3) {
                new_row.data[i] = 1;
            }
            else {
                new_row.data[i] = 0;
            }
        }
    }
    return new_row;
}


// Function to calculate the game of life for the given number of steps
// Function gets the vector of rows for current process, the rank of the process, the size of the world and the number of steps
void calculate(std::vector<row_t>& rows, int rank, int size, int num) {
    row_t up;// the row above the current row
    row_t down;// the row below the current row
    up.data.resize(rows[0].data.size());
    down.data.resize(rows[0].data.size());
    for (int i = 0; i < num; i++) {// cycles over the number of steps
        // first the processes exchanges the rows with their neighbours so they can calculate the new row
        if (rank % 2 == 0) { // even rows send first then wait for resonse while the odd rows wait for the response then send
            // the last process does not have a neighbour below and the first process does not have a neighbour above so
            // they skip a turn
            // since they skip a turn the secont and the second last reciever 2 messages at the same time  to handle
            // this we use the tag so when we are sending the message above we use tag 0 and when we are sending the message
            // below we use tag 1
            if (rank != size - 1) {
                MPI_Send(rows[rows.size() - 1].data.data(), rows[rows.size() - 1].data.size(), MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
                MPI_Recv(down.data.data(), down.data.size(), MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                down.row_number = rows[rows.size() - 1].row_number + 1;
            }
            if (rank != 0) {
                MPI_Send(rows[0].data.data(), rows[0].data.size(), MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
                MPI_Recv(up.data.data(), up.data.size(), MPI_INT, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                up.row_number = rows[0].row_number - 1;
            }
        }
        else {
            if (rank != 0) {
                MPI_Recv(up.data.data(), up.data.size(), MPI_INT, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                up.row_number = rows[0].row_number - 1;
                MPI_Send(rows[0].data.data(), rows[0].data.size(), MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            }
            if (rank != size - 1) {
                MPI_Recv(down.data.data(), down.data.size(), MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                down.row_number = rows[rows.size() - 1].row_number + 1;
                MPI_Send(rows[rows.size() - 1].data.data(), rows[rows.size() - 1].data.size(), MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
            }
        }
        // for easier implementation we add zeroed rows above and below the last and the first row
        if (rank == 0) {
            // the up.data is zeros 
            for (int i = 0; i < rows[0].data.size(); i++) {
                up.data[i] = 0;
            }
            up.row_number = rows[0].row_number - 1;
        }
        if (rank == size - 1) {
            // the down.data is zeros
            for (int i = 0; i < rows[0].data.size(); i++) {
                down.data[i] = 0;
            }
            down.row_number = rows[rows.size() - 1].row_number + 1;
        }
        // the new_rows is the vector of rows of current process after the calculation
        std::vector<row_t> new_rows;
        // if the process has only one row
        // uses only rows up and down to calculate the only row
        if (rows.size() == 1) {
            new_rows.push_back(calc_row(rows[0], up, down));
        }
        //otherwise first is calculated the first row then the middle rows and then the last row
        else {
            new_rows.push_back(calc_row(rows[0], up, rows[1]));
            for (int j = 1; j < rows.size() - 1; j++) { 
                new_rows.push_back(calc_row(rows[j], rows[j - 1], rows[j + 1]));
            }
            new_rows.push_back(calc_row(rows[rows.size() - 1], rows[rows.size() - 2], down));
        }
        rows = new_rows;
    }
}


int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank and size of the world
    int rank, size, numbers_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int len;

    int rows_to_process = 0;
    std::vector<std::vector<int>> numbers;
    std::vector<row_t> my_rows;

    // rank 0 reads the numbers from the file
    if(rank == 0) {
        numbers = read_numbers(argv[1]);

        len = numbers[0].size();
        // broadcasts the number of rows
        numbers_size = numbers.size();
    }

    // Broadcast the number of rows
    MPI_Bcast(&numbers_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast the length of the rows
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process calculates how many rows it will process
    // if the number of rows is not divisible by the number of processes
    // the remaning rows are distributed from start
    rows_to_process = numbers_size / size + ((numbers_size % size > rank) ? 1 : 0);

    // Calculate how many rows each process will get to know how many to send to each 
    int* send_counts = new int[size];
    if (rank == 0)
    {
        for (int i = 0; i < size; i++) {
            send_counts[i] = numbers_size / size + ((numbers_size % size > i) ? 1 : 0);

        }
    }

    // rank 0 sends the rows to the other processes
    if (rank == 0) {
        int i = 0;
        // first it saves the rows it will process
        for (int j = 0; j < rows_to_process; j++) {
            row_t row;
            row.data = numbers[i];
            row.row_number = i;
            my_rows.push_back(row);
            i++;
        }
        // then it sends the rows to the other processes
        for (int j = 1; j < size; j++) {
            for (int k = 0; k < send_counts[j]; k++) {
                MPI_Send(numbers[i].data(), numbers[i].size(), MPI_INT, j, 0, MPI_COMM_WORLD);
                MPI_Send(&i, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                i++;
            }
        }
    }
    else {
        // the other processes recieve the rows from the rank 0
        // and save them to my_rows
        for (int i = 0; i < rows_to_process; i++) {
            int row_number;
            row_t row;
            row.data.resize(len);
            MPI_Recv(row.data.data(), len, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Use data.data() instead of &(data[0])
            MPI_Recv(&row_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            row.row_number = row_number;
            my_rows.push_back(row);
        }
    }

    // the num of steps integer from string
    int num_steps = std::stoi(argv[2]);

    // the process of calculating the game of life
    calculate(my_rows, rank, size, num_steps);

    // barrier to make sure the output is printed in order
    // sometimes does not work :(
    for (int i = 0; i < size; i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == rank)
        {
          print_rows(my_rows, rank);
          fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Finalize MPI
    MPI_Finalize();
}
