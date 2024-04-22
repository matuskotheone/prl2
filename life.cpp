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

struct row {
    std::vector<bool> data;
    int row_number;
};
typedef struct row row_t;

std::vector<std::vector<bool>> read_numbers(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<bool>> numbers;
    while(std::getline(file, line)) {
        std::vector<bool> number;
        for(int i = 0; i < line.size(); i++) {
            number.push_back(line[i] == '1');
        }
        numbers.push_back(number);
    }
    return numbers;
}

void print_rows(std::vector<row_t> rows, int rank) {
    for(int i = 0; i < rows.size(); i++) {
        std::cout << rank << ": " ;
        for(int j = 0; j < rows[i].data.size(); j++) {
            std::cout << (rows[i].data[j] ? "1" : "0");
        }
        std::cout << std::endl;
    }
}

void print_numbers(std::vector<std::vector<bool>> numbers) {
    for(int i = 0; i < numbers.size(); i++) {
        for(int j = 0; j < numbers[i].size(); j++) {
            std::cout << (numbers[i][j] ? "1" : "0");
        }
        std::cout << std::endl;
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
    std::vector<std::vector<bool>> numbers;
    std::vector<row_t> my_rows;

    if(rank == 0) {
        numbers = read_numbers(argv[1]);

        len = numbers[0].size();
        // broadcasts the number of rows
        numbers_size = numbers.size();
    }

    // Broadcast the number of rows
    MPI_Bcast(&numbers_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate how many rows each process will receive
    rows_to_process = numbers_size / size + ((numbers_size % size > rank) ? 1 : 0);

    int* send_counts = new int[size];
    if (rank == 0)
    {
        for (int i = 0; i < size; i++) {
            send_counts[i] = numbers_size / size + ((numbers_size % size > i) ? 1 : 0);

        }
    }

    if (rank == 0) {
        int i = 0;
        for (int j = 0; j < rows_to_process; j++) {
            row_t row;
            row.data = numbers[i];
            row.row_number = i;
            my_rows.push_back(row);
            i++;
        }
        for (int j = 1; j < size; j++) {
            for (int k = 0; k < send_counts[j]; k++) {
                MPI_Send(&numbers[i], numbers[i].size(), MPI_INT, j, 0, MPI_COMM_WORLD);
                MPI_Send(&i, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            }
        }
    }
    else {
        for (int i = 0; i < rows_to_process; i++) {
            int row_number;
            std::vector<bool> data;
            data.resize(len);
            MPI_Recv(&(data[0]), len, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&row_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            row_t row;
            row.data = data;
            row.row_number = row_number;
            my_rows.push_back(row);
        }
    }

    print_rows(my_rows, rank);






    print_rows(my_rows, rank);

    // Finalize MPI
    MPI_Finalize();
}



























