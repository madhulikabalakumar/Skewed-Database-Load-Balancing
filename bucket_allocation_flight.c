#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>

static inline int min(int a, int b) {
    return a < b ? a : b;
}

// Define a struct to store frequency and bucket allocation information
typedef struct {
    int value;          // Value (country or airport code)
    int frequency;      // Frequency of the value
    int startBucket;    // Starting bucket for this value
    int endBucket;      // Ending bucket for this value (exclusive)
} ItemBucket;

// Function prototypes
int round_up(int num, int round);
void init_tuples_and_data_info(uint8_t** data_p, int** dim_ranges_p, int data_info[3], char* filename);

int* calculateFrequencies(uint8_t* data, int tupleNum, int skewColumn, int DIM_RANGE, int NUM_COLS);
int** calculateCombinedFrequencies(uint8_t* data, int tupleNum, int col1, int col2, int* dim_ranges, int NUM_COLS);

ItemBucket* calculateBucketAllocations(int* frequencies, int numValues, int numBuckets, int totalItems);
ItemBucket** calculateCombinedBucketAllocations(int** frequencies, int dim1, int dim2, int tupleNum, int numBuckets);

int round_up(int num, int round) { return (num + round - 1) & (-round); }

void init_tuples_and_data_info(uint8_t** data_p, int** dim_ranges_p, int data_info[3], char* filename) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Could not open file\n");
        return;
    }

    // Count number of columns in CSV
    int num_cols = 1; // Start with 1 to account for the first column
    for (char c = getc(fp); c != EOF; c = getc(fp)) {
        if (c == ',') {
            num_cols++;
        } else if (c == '\n') {
            break;
        }
    }
    data_info[0] = num_cols;

    // Read the dim ranges from the second line of the CSV
    *dim_ranges_p = (int*)malloc(num_cols * sizeof(int));
    int* dim_ranges = *dim_ranges_p;
    fseek(fp, 0, SEEK_SET); // Go back to the beginning of the file
    char line[1024];
    fgets(line, 1024, fp); // Skip the first line (column names)
    fgets(line, 1024, fp); // Read the second line (dim ranges)
    char* token = strtok(line, ",");
    for (int i = 0; i < num_cols; i++) {
        dim_ranges[i] = atoi(token);
        token = strtok(NULL, ",");
    }

    // Count number of lines (i.e., tuples) in CSV file, starting from the third line
    int tuple_num = -2; // Start from -2 to exclude the first two lines
    for (char c = getc(fp); c != EOF; c = getc(fp)) {
        if (c == '\n') {
            tuple_num++;
        }
    }
    data_info[2] = tuple_num;

    // Allocate memory for data
    *data_p = (uint8_t*)malloc(tuple_num * num_cols * sizeof(uint8_t));
    uint8_t* data = *data_p;

    // Read data from CSV file, starting from the third line
    fseek(fp, 0, SEEK_SET); // Go back to the beginning of the file
    fgets(line, 1024, fp); // Skip the first line (column names)
    fgets(line, 1024, fp); // Skip the second line (dim ranges)
    int row = 0;
    while (fgets(line, 1024, fp)) {
        token = strtok(line, ",");
        int col = 0;
        while (token != NULL) {
            data[row * num_cols + col] = (uint8_t)atoi(token);
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(fp);
}

// Function to calculate frequencies of values in the skew column
int* calculateFrequencies(uint8_t* data, int tupleNum, int skewColumn, int DIM_RANGE, int NUM_COLS) {
    int* frequencies = (int*)calloc(DIM_RANGE, sizeof(int));
    for (int i = 0; i < tupleNum; i++) {
        int value = data[i * NUM_COLS + skewColumn];
        frequencies[value]++;
    }
    
    return frequencies;
}

// Function to calculate bucket allocations
ItemBucket* calculateBucketAllocations(int* frequencies, int DIM_RANGE, int tupleNum, int num_buckets) {
    ItemBucket* bucketAllocations = (ItemBucket*)malloc(DIM_RANGE * sizeof(ItemBucket));
    int m = num_buckets; // Start with m = n
    int totalBuckets = 0;

    // Step 1: Initialize all buckets
    for (int i = 0; i < DIM_RANGE; i++) {
        bucketAllocations[i].value = i;
        bucketAllocations[i].frequency = frequencies[i];
        bucketAllocations[i].startBucket = -1; // Initialize to -1
        bucketAllocations[i].endBucket = -1;   // Initialize to -1
    }

    do {
        totalBuckets = 0;
        int lowFreqGroupCount = 0;
        int currentLowFreqBucket = -1;

        // Step 2 & 3: Assign Buckets Based on Frequencies and Group Low Frequency Items
        for (int i = 0; i < DIM_RANGE; i++) {
            if (frequencies[i] == 0) {
                continue; // Skip values with 0 frequency
            }
            if (frequencies[i] * m / tupleNum >= 1) { // Assign buckets to items with frequency * N >= 1
                int bucketsForValue = ceil((double)frequencies[i] * m / tupleNum);
                bucketAllocations[i].startBucket = totalBuckets;
                bucketAllocations[i].endBucket = totalBuckets + bucketsForValue - 1;
                totalBuckets += bucketsForValue;
                lowFreqGroupCount = 0; // Reset low frequency group count
            } else { // Group low-frequency items
                if (lowFreqGroupCount == 0) { // Start a new low-frequency bucket
                    currentLowFreqBucket = totalBuckets;
                    totalBuckets += 1;
                }
                lowFreqGroupCount += frequencies[i];
                bucketAllocations[i].startBucket = currentLowFreqBucket;
                bucketAllocations[i].endBucket = currentLowFreqBucket;
                if (lowFreqGroupCount * m >= tupleNum) {
                    lowFreqGroupCount = 0; // Reset for the next group
                }
            }        
        }

        // Step 4: Adjust m if necessary
        if (totalBuckets > num_buckets) {
            m -= 1; // Decrease m to reduce bucket allocations
        }

    } while (totalBuckets > num_buckets && m > 0);

    // Calculate tuples in each bucket after all allocations are done
    int* tuplesInBucket = (int*)calloc(totalBuckets, sizeof(int));
    for (int i = 0; i < DIM_RANGE; i++) {
        if (frequencies[i] == 0) {
            continue; // Skip values with 0 frequency
        }
        int bucketsAssigned = bucketAllocations[i].endBucket - bucketAllocations[i].startBucket + 1;
        int tuplesPerBucket = frequencies[i] / bucketsAssigned;
        int remainder = frequencies[i] % bucketsAssigned;
        for (int j = bucketAllocations[i].startBucket; j <= bucketAllocations[i].endBucket; j++) {
            tuplesInBucket[j] += tuplesPerBucket + (j - bucketAllocations[i].startBucket < remainder ? 1 : 0);
        }
    }

    // Find the max and min tuples in a bucket
    int maxTuplesInBucket = 0;
    int minTuplesInBucket = INT_MAX;
    int maxBucket = -1; // To track the bucket number with max tuples
    int minBucket = -1; // To track the bucket number with min tuples
    for (int i = 0; i < totalBuckets; i++) {
        if (tuplesInBucket[i] > maxTuplesInBucket) {
            maxTuplesInBucket = tuplesInBucket[i];
            maxBucket = i;
        }
        if (tuplesInBucket[i] < minTuplesInBucket && tuplesInBucket[i] > 0) { // Ignore empty buckets
            minTuplesInBucket = tuplesInBucket[i];
            minBucket = i;
        }
    }

    printf("Max tuples in bucket %d: %d\n", maxBucket, maxTuplesInBucket);
    printf("Min tuples in bucket %d: %d\n", minBucket, minTuplesInBucket);
    printf("Ideal number of tuples per bucket: %d\n", tupleNum / num_buckets);

    free(tuplesInBucket);
    return bucketAllocations;
}

// New function to calculate combined frequencies for tuples from two columns
int** calculateCombinedFrequencies(uint8_t* data, int tupleNum, int col1, int col2, int* dim_ranges, int NUM_COLS) {
    int dim1 = dim_ranges[col1];
    int dim2 = dim_ranges[col2];
    int** frequencies = (int**)malloc(dim1 * sizeof(int*));
    for (int i = 0; i < dim1; ++i) {
        frequencies[i] = (int*)calloc(dim2, sizeof(int));
    }

    for (int i = 0; i < tupleNum; ++i) {
        int val1 = data[i * NUM_COLS + col1];
        int val2 = data[i * NUM_COLS + col2];
        frequencies[val1][val2]++;
    }

    return frequencies;
}

ItemBucket** calculateCombinedBucketAllocations(int** frequencies, int dim1, int dim2, int tupleNum, int numBuckets) {
    // Step 1: Flatten the two-dimensional frequency array
    int totalCombinations = dim1 * dim2; // Correct calculation for total combinations
    int* flattenedFrequencies = (int*)malloc(totalCombinations * sizeof(int));
    int index = 0;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            flattenedFrequencies[index++] = frequencies[i][j];
        }
    }

    // Step 2: Allocate buckets using a more balanced strategy
    ItemBucket* bucketAllocations = calculateBucketAllocations(flattenedFrequencies, totalCombinations, tupleNum, numBuckets);

    // Step 3: Map the allocated buckets back to the two-dimensional allocations array
    ItemBucket** allocations = (ItemBucket**)malloc(dim1 * sizeof(ItemBucket*));
    for (int i = 0; i < dim1; ++i) {
        allocations[i] = (ItemBucket*)malloc(dim2 * sizeof(ItemBucket));
        for (int j = 0; j < dim2; ++j) {
            int flatIndex = i * dim2 + j;
            allocations[i][j].value = j;
            allocations[i][j].frequency = frequencies[i][j];
            allocations[i][j].startBucket = bucketAllocations[flatIndex].startBucket;
            allocations[i][j].endBucket = bucketAllocations[flatIndex].endBucket;
        }
    }

    // Clean up
    free(flattenedFrequencies);
    free(bucketAllocations);

    return allocations;
}


int main(int argc, char* argv[]) {
    // Initialize tuple data
    uint8_t* data;
    int* dim_ranges;
    int data_info[3];
    init_tuples_and_data_info(&data, &dim_ranges, data_info, argv[1]);
    const int NUM_COLS = data_info[0];
    const int tuple_num = data_info[2];

    int num_buckets = 128;

    FILE *fp = fopen("bucket_allocations_flight_discrete.csv", "w");
    fprintf(fp, "Column,Value,Frequency,StartBucket,EndBucket\n");

    // Loop through each column to calculate frequencies and bucket allocations
    for (int col = 0; col < NUM_COLS; col++) {
        int numValues = dim_ranges[col];
        int* frequencies = calculateFrequencies(data, tuple_num, col, numValues, NUM_COLS);
        ItemBucket* bucketAllocations = calculateBucketAllocations(frequencies, numValues, tuple_num, num_buckets);

        // Save bucket allocations to CSV
        for (int i = 0; i < numValues; i++) {
            printf("%d,%d,%d,%d,%d\n", col, i, frequencies[i],
                    bucketAllocations[i].startBucket, bucketAllocations[i].endBucket);
            fprintf(fp, "%d,%d,%d,%d,%d\n", col, i, frequencies[i],
                    bucketAllocations[i].startBucket, bucketAllocations[i].endBucket);
        }

        // Clean up
        free(frequencies);
        free(bucketAllocations);
    }

    fclose(fp);

    FILE *fp2 = fopen("bucket_allocations_flight_discrete_grouped.csv", "w");
    if (fp2 == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    fprintf(fp2, "Column1,Column2,Value1,Value2,Frequency,StartBucket,EndBucket\n");

    // Process combinations of two columns
    for (int col1 = 0; col1 < NUM_COLS; col1++) {
        for (int col2 = col1 + 1; col2 < NUM_COLS; col2++) {
            // Calculate combined frequencies for each tuple
            int** combinedFrequencies = calculateCombinedFrequencies(data, tuple_num, col1, col2, dim_ranges, NUM_COLS);

            // Calculate combined bucket allocations for the column pair
            ItemBucket** combinedAllocations = calculateCombinedBucketAllocations(combinedFrequencies, dim_ranges[col1], dim_ranges[col2], tuple_num, num_buckets);

            // Save combined bucket allocations to CSV
            for (int i = 0; i < dim_ranges[col1]; i++) {
                for (int j = 0; j < dim_ranges[col2]; j++) {
                    if (combinedAllocations[i][j].frequency > 0) { // Skip zero frequencies
                        fprintf(fp2, "%d,%d,%d,%d,%d,%d,%d\n", col1, col2, i, j,
                                combinedAllocations[i][j].frequency,
                                combinedAllocations[i][j].startBucket, 
                                combinedAllocations[i][j].endBucket);
                        printf("%d,%d,%d,%d,%d,%d,%d\n", col1, col2, i, j,
                                combinedAllocations[i][j].frequency,
                                combinedAllocations[i][j].startBucket, 
                                combinedAllocations[i][j].endBucket);
                    }
                }
            }

            // Free combined allocations and frequencies
            for (int i = 0; i < dim_ranges[col1]; i++) {
                free(combinedAllocations[i]);
                free(combinedFrequencies[i]);
            }
            free(combinedAllocations);
            free(combinedFrequencies);
        }
    }

    // Clean up
    fclose(fp2);
    free(data);
    free(dim_ranges);

    return 0;
}

