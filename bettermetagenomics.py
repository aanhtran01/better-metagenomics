# -*- coding: utf-8 -*-
'''
Goal: perform metagenomic analysis on a sample.

Input: genome sequences which maybe in a sample and reads from the sample, certain number of genomes are present in each sample
Output: headers of the reads along with the source genome of the read. predicted sources are compared to the true sources and the score is based on how many reads are classified correctly.
'''

import argparse
from copy import deepcopy
import os
import re
import zipfile
import random
import itertools

def extract_reference_genome(reference_file):
    # Open the reference file and extract the DNA sequence
    with open(reference_file) as f:
        # Read the file
        lines = f.readlines()
        # Remove newlines and any leading/trailing spaces
        lines = [line.strip() for line in lines]
        # Concatenate the DNA sequence lines
        reference_genome = "".join(lines[1:])

    return reference_genome

def extract_genome_name(file_path):
    with open(file_path, 'r') as file:
        header_line = file.readline().strip()
        match = re.search(r">Genome_Number_(\d+)", header_line)
        if match:
            genome_id = match.group(1)
            return f"Genome_Number_{genome_id}"
        else:
            return None
'''
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        donor_id = None
        donor_seq = ""
        for line in f:
            if line.startswith(">"):
                if donor_id is not None:
                    yield (donor_id, donor_seq)
                donor_id = line.strip()[1:]
                donor_seq = ""
            else:
                donor_seq += line.strip()
        if donor_id is not None:
            yield (donor_id, donor_seq)
'''

# Function to read in donor reads fasta file
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        donor_id = None
        donor_seq = ""
        read_count = 0  # Counter for the number of reads
        for line in f:
            if line.startswith(">"):
                if donor_id is not None:
                    yield (donor_id, donor_seq)
                    read_count += 1
                    if read_count >= 100000:
                        break  # Stop after 100,000 reads
                donor_id = line.strip()[1:]
                donor_seq = ""
            else:
                donor_seq += line.strip()
        if donor_id is not None and read_count < 100000:
            yield (donor_id, donor_seq)
            

def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])

#function to align a read to the reference genome
def align_read_to_genome(read, reference_genome, k, reference_kmers):

    # Step 1: Break the read into k-mers
    kmers = [read[i:i+k] for i in range(len(read)-k+1)]
    #first_kmer = kmers[0]


    # Step 2: Search for matches in the BurrowsWheeler index and extend the alignment
    best_match = None
    best_score = float('inf')
    for i, kmer in enumerate(kmers):
        positions = reference_kmers.get(kmer)

        if positions is not None:  # Check if positions exist
            for pos in positions:
                # extend the alignment
                offset = i
                alignment_start = pos - offset
                alignment_end = alignment_start + len(read)
                if alignment_start < 0 or alignment_end > len(reference_genome):
                    continue  # alignment out of bounds, skip to next position
                ref_sequence = reference_genome[alignment_start:alignment_end]
                score = HammingDistance(read, ref_sequence)

                # check if this is the best match so far
                if score < best_score:
                    best_score = score
                    best_match = alignment_start

    return best_match, best_score

def align_all_reads_to_genome(donor_reads, reference_genome, reference_kmers, k, genome_name, existing_results=None):
    results = existing_results if existing_results is not None else []

    for read_id, read_seq in donor_reads:
        best_match, best_score = align_read_to_genome(read_seq, reference_genome, k, reference_kmers)

        # Check if read_id already exists in results
        existing_entry = next((entry for entry in results if entry['donor_read_id'] == read_id), None)

        if existing_entry is None or best_score < existing_entry['best_score']:
            result = {
                'donor_read_id': read_id,
                'sequence': read_seq,
                'best_match': best_match,
                'best_score': best_score,
                'genome_id': genome_name
            }

            if existing_entry is None:
                results.append(result)
            else:
                existing_entry.update(result)

    return results

def pseudo_align(read_id, read, reference_minimizers, read_minimizers):
    # Step 1: Retrieve the read minimizers from the minimizer dictionary
    read_minimizers = read_minimizers.get(read_id, [])
    count = 0

    for i, minimizer in enumerate(read_minimizers):
        if minimizer in reference_minimizers:
            count += 1  # Matching minimizer found

    return count


def align_sampled_reads_to_genome(donor_reads, reference_minimizers, minimizer_size, window_size, genome_name, read_minimizers, existing_scores=None):
    scores = existing_scores if existing_scores is not None else []

    for genome_id in genome_name:
        count = 0
        for read_id, read_seq in donor_reads:
            count += pseudo_align(read_id, read_seq, reference_minimizers, read_minimizers)
        if count > 0:
            scores.append({'count': count, 'genome_id': genome_name})
        unique_genomes = set(score['genome_id'] for score in scores)
    top_scores = []

    for genome in unique_genomes:
        max_count = max(score['count'] for score in scores if score['genome_id'] == genome)
        top_scores.append({'count': max_count, 'genome_id': genome})

    return top_scores

# Create an argument parser
parser = argparse.ArgumentParser(description='Extract files from a zip archive.')

# Add the zip file path argument
parser.add_argument('zip_file', type=str, help='Path to the zip file.')

# Parse the command-line arguments
args = parser.parse_args()

# Get the current working directory
folder_path = os.getcwd()

# Extract the files from the zip archive to the current working directory
with zipfile.ZipFile(args.zip_file, 'r') as zip_ref:
    zip_ref.extractall(folder_path)
    
# Get all the files in the extracted folder
file_list = os.listdir(folder_path)

# Filter the files based on the condition
genome_list = [file for file in file_list if "genome" in file]
#print(genome_list)

# Filter the files based on the condition
read_list = [file for file in file_list if "reads.fasta" in file]

# Extract the first file from the list
read_file = read_list[0]

def generate_read_minimizers(read, read_id, window_size, minimizer_size, read_minimizers):
    read_minimizers[read_id] = set()
    for i in range(len(read) - window_size + 1):
        window = read[i:i + window_size]
        min_substring = min(window[i:i + minimizer_size] for i in range(window_size - minimizer_size + 1))
        read_minimizers[read_id].add(min_substring)

def generate_minimizers(sequence, window_size, minimizer_size):
    minimizer_dict = {}

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        min_substring = min(window[i:i + minimizer_size] for i in range(window_size - minimizer_size + 1))

        if min_substring in minimizer_dict:
            minimizer_dict[min_substring].append(i)
        else:
            minimizer_dict[min_substring] = [i]

    return minimizer_dict

def generate_reference_kmers(reference_genome, kmer_size):
    reference_kmers = {}
    for i in range(len(reference_genome) - kmer_size + 1):
        kmer = reference_genome[i:i+kmer_size]
        position = i  # Adjust the position to be 1-indexed
        if kmer in reference_kmers:
            reference_kmers[kmer].append(position)
        else:
            reference_kmers[kmer] = [position]
    return reference_kmers

window_size = 21
minimizer_size = 15

read_path = os.path.join(folder_path, read_file)
donor_reads = list(read_fasta_file(read_path))  # Convert generator to a list

read_minimizers = {}
for read_id, read_seq in donor_reads:
    generate_read_minimizers(read_seq, read_id, window_size, minimizer_size, read_minimizers)

sample_size = int(200)

# Randomly sample reads from donor_reads
sampled_reads = random.sample(donor_reads, sample_size)


sample_results = []

# Process each file one by one
for genome_name in genome_list:
    file_path = os.path.join(folder_path, genome_name)

    # Perform your processing on each file
    # Call the function to extract the DNA sequence
    reference_genome = extract_reference_genome(file_path)

    reference_minimizers = generate_minimizers(reference_genome, window_size, minimizer_size)
    #genome_id = extract_genome_name(file_path)

    # Use the sampled_reads for further processing
    sample_alignment = align_sampled_reads_to_genome(sampled_reads, reference_minimizers , minimizer_size, window_size, genome_name, read_minimizers, sample_results)
    #print(sample_alignment)
    sample_results = sample_alignment
    #sample_results.extend(sample_alignment)  # Append the new results to the existing ones

k=25
filtered_results = [item for item in sample_results if item['count'] > 18]
#filtered_results = [item for item in sample_results if item['count'] > 70]
results = []
# Process each file one by one
for genome_id in filtered_results:
    genome_id_value = genome_id['genome_id']  # Access the value of 'genome_id'
    path = os.path.join(folder_path, genome_id_value)

    # Perform your processing on each file
    # Call the function to extract the DNA sequence
    reference_genome = extract_reference_genome(path)

    reference_kmers = generate_reference_kmers(reference_genome, k)
    genome_id = extract_genome_name(path)

    #read_path = os.path.join(folder_path, read_file)
    #donor_reads = read_fasta_file(read_path)

    # Pass the existing results dictionary to the alignment function
    alignment = align_all_reads_to_genome(donor_reads, reference_genome, reference_kmers, k, genome_id, results)

    results = alignment
output_file = "predictions.txt"  # Specify the path to the output file

# Write results to the output file
with open(output_file, 'w') as file:
    for result in results:
        file.write(">")
        line = f"{result['donor_read_id']}\t{result['genome_id']}\n"
        file.write(line)
