#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <Rcpp.h>
#include <limits.h>
#include <unistd.h>

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
int FastaSlice(std::string query, std::string target, int slicesize, int leaps) {
    char seq_id[200], seq_str[101], ch, command[100], db_seqs[50];
    int i, seq_id_flag = 0, seq_id_len = 0, seq_str_len = 0, success = 0,
            size = 0, unk = 0, newline = 0, skips = 0;
    long int line_num = 0, seq_pos = 0;
    const char *query_genome = query.c_str(), *target_database = target.c_str();
    vector<int> qchrom_lengths, dbchrom_lengths;
    vector <string> db_seq_names, q_seq_names;
    FILE *query_ptr, *queryout_ptr, *results_ptr, *db_ptr, *chrom_ptr;
    time_t begin = time(NULL);

    query_ptr = fopen(query_genome, "r");
    queryout_ptr = fopen("query_v02.fsa", "w");
    fclose(queryout_ptr);
    results_ptr = fopen("results_v03.out", "w");
    fclose(results_ptr);

    printf("Shaking");

    while ((ch = getc(query_ptr)) != EOF)//loop through the fasta file
    {
        if (ch == '\n') //count the number of lines in the file
        {
            line_num++;
        }

        if (ch == '>')//detect beginning of sequence name
        {
            if (seq_pos != 0)
                qchrom_lengths.push_back(seq_pos);

            seq_pos = 0;
            seq_id_flag = 1;//set flags
            seq_id_len = 1;
            ch = getc(query_ptr); //exclude the carrot
            while (seq_id_flag == 1 && ch != EOF)//read in sequence name string
            {
                seq_id[seq_id_len - 1] = ch;
                seq_id_len++;
                ch = getc(query_ptr);
                if (ch == '\n') {
                    seq_id_flag = 0;//remove flags
                    seq_id[seq_id_len] = '\0';
                    line_num++;
                    q_seq_names.push_back(seq_id);
                }
            }
        }

        if (ch != '>' && ch != '\n' && seq_id_flag == 0) //if no longer a sequence name string then start reading bases
        {
            while (skips < leaps && ch != '>')
            {
                if ((ch = getc(query_ptr)) != '\n')
                {
                    seq_pos++;
                    skips++;
                }
                if (ch == '\n')
                {
                    line_num++;
                }
                if (ch == '>')
                {
                    fseek(query_ptr, -1, SEEK_CUR);
                }
            }
            while (seq_str_len < slicesize && ch != '>' && ch != EOF && skips >= leaps)//split into 100 base blocks
            {
                if ((ch = getc(query_ptr)) != '\n' && ch != '>' && ch != EOF)
                {
                    seq_pos++;
                    seq_str_len++;
                    seq_str[seq_str_len - 1] = ch;
                }
                if (ch == '\n')
                    line_num++;
                if (ch == 'N')
                    unk++;
                if (ch == '>')
                    fseek(query_ptr, -1, SEEK_CUR);
            }

            seq_str[seq_str_len] = '\0';

            if (unk < slicesize / 2 && strlen(seq_str) > 60)
            {
                queryout_ptr = fopen("query_v02.fsa", "a");
                fprintf(queryout_ptr, ">%li_%s\n", seq_pos, seq_id);
                fprintf(queryout_ptr, "%s\n", seq_str);
                fclose(queryout_ptr);
            }
        }
        skips = 0;
        seq_str_len = 0;
        unk = 0;
        if (!(line_num % 1000))
            printf(".");
    }

    qchrom_lengths.push_back(seq_pos);

    sprintf(command,
            "blastn -db %s -query query_v02.fsa  -outfmt \"6 qseqid sstart sseqid length pident score evalue\" -max_target_seqs 1 -out results_v03.out",
            target_database);
    success = system(command);

    printf("\nBaking");
    results_ptr = fopen("results_v03.out", "r+");
    newline = 1;
    line_num = 0;
    while ((ch = getc(results_ptr)) != EOF) {
        if (newline && ch == '_') {
            fseek(results_ptr, -1, SEEK_CUR);
            fprintf(results_ptr, "\t");
            newline = 0;
        }
        if (ch == '\n') {
            newline = 1;
            line_num++;
            if (!(line_num % 500))
                printf(".");
        }
    }
    fclose(results_ptr);
    fclose(query_ptr);

    //pt 2
    db_ptr = fopen(target_database, "r");
    seq_pos = 0;
    while ((ch = getc(db_ptr)) != EOF) {
        if (ch == '>')//detect beginning of sequence name
        {
            if (seq_pos != 0)
                dbchrom_lengths.push_back(seq_pos);
            seq_pos = 0;
            success = fscanf(db_ptr, "%s", db_seqs);
            db_seq_names.push_back(db_seqs);
            while ((ch = getc(db_ptr)) != '\n') {
            }
        }
        seq_pos++;
    }
    dbchrom_lengths.push_back(seq_pos);
    fclose(db_ptr);

    chrom_ptr = fopen("chromlength.out", "w");

    for (i = 0; i < qchrom_lengths.size(); i++) {
        fprintf(chrom_ptr, "%s\t", q_seq_names[i].c_str());
    }
    fprintf(chrom_ptr, "\n");
    for (i = 0; i < qchrom_lengths.size(); i++) {
        fprintf(chrom_ptr, "%i\t", qchrom_lengths[i]);
    }
    fprintf(chrom_ptr, "\n");
    for (i = 0; i < dbchrom_lengths.size(); i++) {
        fprintf(chrom_ptr, "%i\t", dbchrom_lengths[i]);
    }
    fprintf(chrom_ptr, "\n");
    for (i = 0; i < db_seq_names.size(); i++) {
        fprintf(chrom_ptr, "%s\t", db_seq_names[i].c_str());
    }

    fclose(chrom_ptr);

    time_t end = time(NULL);
    double time_spent = (double) (end - begin) / 60;
    printf("\nComplete!\nTime: %lf mins", time_spent);
    return 0;
}
