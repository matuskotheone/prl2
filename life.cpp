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
    std::vector<int> data;
    int row_number;
};
typedef struct row row_t;

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

void print_row(row_t row) {
    for(int j = 0; j < row.data.size(); j++) {
        std::cout << (row.data[j] ? "1" : "0");
    }
    std::cout << std::endl;
}

void print_rows(std::vector<row_t> rows, int rank) {
    for(int i = 0; i < rows.size(); i++) {
        std::cout << rank << ": " ;
        print_row(rows[i]);
    }
}

void print_numbers(std::vector<std::vector<int>> numbers) {
    for(int i = 0; i < numbers.size(); i++) {
        for(int j = 0; j < numbers[i].size(); j++) {
            std::cout << (numbers[i][j] ? "1" : "0");
        }
        std::cout << std::endl;
    }
}

row_t calc_row(row_t curr, row_t up, row_t down) {
    row_t new_row;
    new_row.data.resize(curr.data.size());
    new_row.row_number = curr.row_number;
    for (int i = 0; i < curr.data.size(); i++) {
        int count = 0;
        if (i == 0) {
            count += up.data[i];
            count += up.data[i + 1];
            count += curr.data[i + 1];
            count += down.data[i];
            count += down.data[i + 1];
        }
        else if (i == curr.data.size() - 1) {
            count += up.data[i - 1];
            count += up.data[i];
            count += curr.data[i - 1];
            count += down.data[i - 1];
            count += down.data[i];
        }
        else {
            count += up.data[i - 1];
            count += up.data[i];
            count += up.data[i + 1];
            count += curr.data[i - 1];
            count += curr.data[i + 1];
            count += down.data[i - 1];
            count += down.data[i];
            count += down.data[i + 1];
        }
        if (curr.data[i] == 1) {
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
        else {
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


void calculate(std::vector<row_t>& rows, int rank, int size, int num) {
    row_t up;
    row_t down;
    up.data.resize(rows[0].data.size());
    down.data.resize(rows[0].data.size());
    for (int i = 0; i < num; i++) {
        if (rank % 2 == 0) {
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
        std::vector<row_t> new_rows;
        if (rows.size() == 1) {
            new_rows.push_back(calc_row(rows[0], up, down));
        }
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
                MPI_Send(numbers[i].data(), numbers[i].size(), MPI_INT, j, 0, MPI_COMM_WORLD);
                MPI_Send(&i, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                i++;
            }
        }
    }
    else {
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

    int num_steps = std::stoi(argv[2]);

    calculate(my_rows, rank, size, num_steps);
    int i = 0;

    for (int i = 0; i < size; i++)
    {
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



























