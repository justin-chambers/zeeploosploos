#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <exception>

#define MAX_PROCESSORS  8
#define MIN_PROBLEM_SIZE    120
#define MAX_PROBLEM_SIZE    10000
#define MIN_SAMPLE_SIZE     60
#define MAX_SAMPLE_SIZE     120
#define GENERATE_SEQ    0
#define GENERATE_PARA   1
#define GENERATE_FULL_SET   0
#define GENERATE_SAMPLE_SET 1

void write_slurm_shell(const char * shell_name,
                       const char * bin_name,
                       const int problem_size,
                       const int num_processors,
                       const char * sample_name,
                       const char * log_name);

void write_slurm_shell_with_files(const char * shell_name,
                                  const char * bin_name,
                                  const int num_processors,
                                  const char * input_one,
                                  const char * input_two,
                                  const char * output,
                                  const char * log_name);

int main() {
    std::vector<std::vector<int>> problem_sets;
    std::vector<int> set_1;

#if GENERATE_FULL_SET
    for(int i=MIN_PROBLEM_SIZE; i<=MAX_PROBLEM_SIZE; i=i+60) {
        set_1.push_back(i);
    }
#endif

#if GENERATE_SAMPLE_SET
    for(int i=MIN_SAMPLE_SIZE; i<=MAX_SAMPLE_SIZE; i=i+60) {
        set_1.push_back(i);
    }
#endif

#if GENERATE_SEQ
    std::vector<int> processors = { 1 };
    std::string bin = "../bin/mm_seq";
    std::string sample = "../results/s";
    std::string log = "log_s.csv";
    std::string shell = "s";
#endif

#if GENERATE_PARA
    std::vector<int> processors = { 4, 9, 16, 25, 36 };
    std::string bin = "../bin/mm_para";
    std::string sample = "../results/p";
    std::string log = "log_p.csv";
    std::string shell = "p";
#endif
    problem_sets.push_back(set_1);
    for(int k=0; k<(int)problem_sets.size(); ++k) {
        std::vector<int> set = problem_sets[k];
        for(int i=0; i<(int)set.size();++i) {
            for(int j=0; j<(int)processors.size();++j) {
                std::string file1 = "../tests/a";
                std::string file2 = "../tests/b";
                file1 = file1 + std::to_string(set[i]) + ".txt";
                file2 = file2 + std::to_string(set[i]) + ".txt";
                std::string output = sample;
                output = output + "_" + std::to_string(set[i]) + "_" + std::to_string(processors[j]) + ".txt";
                shell = shell + std::to_string(set[i]) + "_" + std::to_string(processors[j]) + ".sh";
                write_slurm_shell_with_files(shell.c_str(),
                                             bin.c_str(),
                                             processors[j],
                                             file1.c_str(),
                                             file2.c_str(),
                                             output.c_str(),
                                             log.c_str());
#if GENERATE_SEQ
                shell = "s";
#endif
#if GENERATE_PARA
                shell = "p";
#endif
            }
        }
        set.clear();
    }
    return 0;
}

void write_slurm_shell(const char * shell_name,
                       const char * bin_name,
                       const int problem_size,
                       const int num_processors,
                       const char * sample_name,
                       const char * log_name) {
    std::fstream fout;
    try {
        fout.open(shell_name, std::ios::out);
    } catch ( std::exception e ) {
        std::cout << e.what() << std::endl;
    }

    // write SBATCH header...
    fout << "#!/bin/bash\n\n";

    #if GENERATE_SEQ
    fout << "#SBATCH -n 1\n";
    fout << "#SBATCH -N 1\n";
#endif

#if GENERATE_PARA
    fout << "#SBATCH -n " << num_processors <<"\n";
    fout << "#SBATCH -N " << (num_processors/MAX_PROCESSORS)+1 << "\n";
#endif

    fout << "#SBATCH --cpus-per-task=1\n";
    fout << "#SBATCH --mem=2048MB\n";
    fout << "#SBATCH --time=00:10:00\n";
    fout << "#SBATCH --mail-user=justinchambers@nevada.unr.edu\n";
    fout << "#SBATCH --mail-type=ALL\n\n";

    // write SBATCH command...
    fout << "srun "
         << bin_name << " "
         << problem_size << " "
         << num_processors << " "
         << sample_name << problem_size << "_" << num_processors << ".out"
         << " >> " << log_name << std::endl;
    fout.close();
}

void write_slurm_shell_with_files(const char * shell_name,
                                  const char * bin_name,
                                  const int num_processors,
                                  const char * input_one,
                                  const char * input_two,
                                  const char * output,
                                  const char * log_name) {
    std::fstream fout;
    try {
        fout.open(shell_name, std::ios::out);
    } catch ( std::exception e ) {
        std::cout << e.what() << std::endl;
    }

    // write SBATCH header...
    fout << "#!/bin/bash\n\n";

#if GENERATE_SEQ
    fout << "#SBATCH -n 1\n";
    fout << "#SBATCH -N 1\n";
#endif

#if GENERATE_PARA
    fout << "#SBATCH -n " << num_processors <<"\n";
    fout << "#SBATCH -N " << (num_processors/MAX_PROCESSORS)+1 << "\n";
#endif

    fout << "#SBATCH --cpus-per-task=1\n";
    fout << "#SBATCH --mem=2048MB\n";
    fout << "#SBATCH --time=00:10:00\n";
    fout << "#SBATCH --mail-user=justinchambers@nevada.unr.edu\n";
    fout << "#SBATCH --mail-type=ALL\n\n";

    // write SBATCH command...
    fout << "srun "
         << bin_name << " "
         << num_processors << " "
         << input_one << " "
         << input_two << " "
         << output << " >> " << log_name << std::endl;
    fout.close();
}