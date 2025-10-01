import os
from collections import defaultdict

def parse_maf_and_merge_with_N(maf_path, spacer_length=100):
    sequences = defaultdict(list)
    spacer = 'N' * spacer_length

    # Derive output FASTA path
    base_name = os.path.splitext(os.path.basename(maf_path))[0]
    output_fasta = os.path.join(os.path.dirname(maf_path), f"{base_name}.fasta")

    # First pass: find maximum mult value
    max_mult = 0
    with open(maf_path, 'r') as infile:
        for line in infile:
            if line.startswith('a ') and 'mult=' in line:
                mult_val = int(line.strip().split('mult=')[-1])
                max_mult = max(max_mult, mult_val)

    # Second pass: extract and concatenate only full alignments
    with open(maf_path, 'r') as infile:
        block_lines = []
        include_block = False

        for line in infile:
            if line.startswith('a '):
                # Process previous block
                if block_lines and include_block:
                    for bline in block_lines:
                        parts = bline.strip().split()
                        seq_id = parts[1]
                        aligned_seq = parts[6]
                        sequences[seq_id].append(aligned_seq)
                    # Add spacer
                    for seq_id in sequences:
                        sequences[seq_id].append(spacer)

                # Start new block
                block_lines = []
                include_block = f"mult={max_mult}" in line

            elif line.startswith('s ') and include_block:
                block_lines.append(line)

        # Final block after file ends
        if block_lines and include_block:
            for bline in block_lines:
                parts = bline.strip().split()
                seq_id = parts[1]
                aligned_seq = parts[6]
                sequences[seq_id].append(aligned_seq)
            for seq_id in sequences:
                sequences[seq_id].append(spacer)

    # Finalize sequences
    final_sequences = {}
    for seq_id, seq_list in sequences.items():
        if seq_list and seq_list[-1] == spacer:
            seq_list.pop()
        final_sequences[seq_id] = ''.join(seq_list)

    # Write to FASTA
    with open(output_fasta, 'w') as outfile:
        for seq_id, full_seq in final_sequences.items():
            outfile.write(f">{seq_id}\n")
            for i in range(0, len(full_seq), 80):
                outfile.write(full_seq[i:i+80] + '\n')

    print(f"Output written to: {output_fasta}")

#Example usage:
parse_maf_and_merge_with_N(r"C:\Users\User\Documents\Bioinformatics_Year3_Sem1\FYP\mugsy_output\s_rhizophila\tmp.maf")


# from collections import defaultdict
# from math import ceil

# def parse_maf_with_partial_blocks(maf_file, output_fasta, spacer_length=100, min_genome_fraction=0.9):
    # sequences = defaultdict(list)
    # spacer = 'N' * spacer_length
    # all_genomes = set()

    # # First pass: collect all genome IDs
    # with open(maf_file, 'r') as infile:
        # for line in infile:
            # if line.startswith('s '):
                # seq_id = line.strip().split()[1]
                # all_genomes.add(seq_id)

    # total_genomes = len(all_genomes)
    # print(f"Total genomes found: {total_genomes}")

    # # Second pass: process blocks
    # with open(maf_file, 'r') as infile:
        # block_lines = []
        # current_block_genomes = set()

        # for line in infile:
            # if line.startswith('a '):
                # # Process the previous block
                # if block_lines:
                    # if len(current_block_genomes) / total_genomes >= min_genome_fraction:
                        # block_seqs = {}
                        # block_length = None
                        # for bline in block_lines:
                            # parts = bline.strip().split()
                            # seq_id = parts[1]
                            # aligned_seq = parts[6]
                            # block_seqs[seq_id] = aligned_seq
                            # if block_length is None:
                                # block_length = len(aligned_seq)
                        # # Add aligned or padded gaps to each genome
                        # for genome in all_genomes:
                            # if genome in block_seqs:
                                # sequences[genome].append(block_seqs[genome])
                            # else:
                                # sequences[genome].append('-' * block_length)
                        # # Add spacer
                        # for genome in all_genomes:
                            # sequences[genome].append(spacer)

                # # Start new block
                # block_lines = []
                # current_block_genomes = set()

            # elif line.startswith('s '):
                # block_lines.append(line)
                # seq_id = line.strip().split()[1]
                # current_block_genomes.add(seq_id)

        # # Final block after file ends
        # if block_lines:
            # if len(current_block_genomes) / total_genomes >= min_genome_fraction:
                # block_seqs = {}
                # block_length = None
                # for bline in block_lines:
                    # parts = bline.strip().split()
                    # seq_id = parts[1]
                    # aligned_seq = parts[6]
                    # block_seqs[seq_id] = aligned_seq
                    # if block_length is None:
                        # block_length = len(aligned_seq)
                # for genome in all_genomes:
                    # if genome in block_seqs:
                        # sequences[genome].append(block_seqs[genome])
                    # else:
                        # sequences[genome].append('-' * block_length)
                # for genome in all_genomes:
                    # sequences[genome].append(spacer)

    # # Finalize sequences
    # final_sequences = {}
    # for seq_id, seq_list in sequences.items():
        # if seq_list and seq_list[-1] == spacer:
            # seq_list.pop()
        # final_sequences[seq_id] = ''.join(seq_list)

    # # Write to FASTA
    # with open(output_fasta, 'w') as outfile:
        # for seq_id, full_seq in final_sequences.items():
            # outfile.write(f">{seq_id}\n")
            # for i in range(0, len(full_seq), 80):
                # outfile.write(full_seq[i:i+80] + '\n')

# # Usage
# parse_maf_with_partial_blocks("mga_non_pathogen.maf", "aligned_non_pathogen_mugsy_partial90.fasta")

