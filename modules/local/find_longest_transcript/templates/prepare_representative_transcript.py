#!/usr/bin/env python3

"""
Calculates a table of CDS and exon lengths for each transcript and validates
or autoselects representative transcripts.

Validates user-provided transcripts or automatically selection of the longest
transcript per gene based on CDS and exon length.

Outputs:
- A list of selected transcript IDs (.txt), if not provided
- A transcript length index file (.fai) with transcript ID and length
- A transcript GTF
"""

import platform
import argparse
from sys import exit
import pyranges as pr
import pandas as pd
import warnings
import csv
import logging
import os


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")

def parse_gtf(gtf_path):
    """
    Reads a GTF file and checks that every gene has at least one transcript.

    Returns:
    - df_gtf: full parsed GTF as a DataFrame
    """
    df_gtf = pr.read_gtf(gtf_path, as_df=True)
    input_genes = set(df_gtf.gene_id)

    # Identify genes that have at least one associated transcript
    transcript_rows = df_gtf[df_gtf['Feature'] == 'transcript']
    genes_with_transcripts = set(transcript_rows.gene_id)

    # Check for genes with no transcripts
    missing_genes = input_genes - genes_with_transcripts
    if missing_genes:
        top = sorted(missing_genes)[:10]
        more = f" and {len(missing_genes) - 10} more" if len(missing_genes) > 10 else ""
        logging.error(f"Some genes have no transcript entries in the GTF: {top}{more}")
        raise ValueError(
            "ERROR: Some genes in the GTF have no associated transcript entries. "
            "Please ensure the GTF includes transcript-level features for each gene."
        )

    return df_gtf


def parse_gtf_and_calculate_lengths(gtf_path):
    """
    Reads a GTF file and calculates CDS and exon lengths per transcript.

    Returns:
    - df_gtf: full parsed GTF as a DataFrame (with added 'length' column)
    - df_txlengths: DataFrame with gene_id, transcript_id, cds_length, exon_length
    - input_genes: set of all gene_ids present in the GTF
    """
    df_gtf = parse_gtf(gtf_path)
    df_gtf['length'] = df_gtf['End'] - df_gtf['Start']
    input_genes = set(df_gtf.gene_id)

    # Calculate CDS and exon lengths per transcript
    cds_sums = df_gtf[df_gtf['Feature'] == 'CDS'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    exon_sums = df_gtf[df_gtf['Feature'] == 'exon'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    tx_sums = df_gtf[df_gtf['Feature'] == 'transcript'].copy().set_index(['gene_id', 'transcript_id'])['length']
    df_txlengths = pd.concat([cds_sums, exon_sums, tx_sums], axis=1, keys=['cds_length', 'exon_length', 'unspliced_length']).reset_index()

    # Warn if any transcript is missing exon length
    missing_exons = df_txlengths[df_txlengths.exon_length.isna()]
    if not missing_exons.empty:
        warnings.warn("Some transcripts do not have exon length.", RuntimeWarning)
        print("Transcripts missing exon length:")
        print(missing_exons[['gene_id', 'transcript_id']])

    # Fill missing values with 0
    df_txlengths[['cds_length', 'exon_length', 'unspliced_length']] = df_txlengths[['cds_length', 'exon_length', 'unspliced_length']].fillna(0)

    # Check if all genes have at least one transcript
    missing_genes = input_genes - set(df_txlengths.gene_id)
    if len(missing_genes) > 0:
        logging.error("Some genes don't have transcripts associated with them.")
        raise ValueError(
            "ERROR: Some genes have no associated transcripts. The GTF may be corrupted.\n"
            "Please provide a valid GTF with at least one transcript per gene, or run with --filter_gtf_by_transcripts false."
        )


    return df_gtf, df_txlengths, input_genes


def main(process_name, gtf, transcript, output):
    # Dump versions file
    dump_versions(process_name)

    # Logging
    log_file = f"{output}.log"
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load the genome annotation file and calculate CDS and exon lengths per transcript
    df_gtf, df_txlengths, input_genes = parse_gtf_and_calculate_lengths(gtf)
    df_txlengths.to_csv(f"{os.path.basename(gtf)}.lengths.tsv", index=None, sep='\t', quoting=csv.QUOTE_NONE)
    logging.info(f"Total genes in genome GTF: {len(input_genes)}")

    # If user specified transcripts are provided, use these as representative transcripts; else find longest CDS / exon transcript.
    if transcript != "":
        logging.info(f"Using the transcript IDs in: {transcript}")
        with open(transcript, "r") as file:
            transcript_ids = [line.strip() for line in file]
            # Remove potential empty strings caused by trailing newline
            transcript_ids = [x for x in transcript_ids if x]
            transcript_ids_sorted = transcript_ids

        # Check if the transcript file was empty or only contained blank lines
        if len(transcript_ids) == 0:
            logging.error("Transcript list file is empty.")
            raise ValueError(
                "ERROR: The transcript list provided is empty. "
                "Please provide a non-empty file with transcript IDs."
            )

        # Transcript list validation: check if all user-provided IDs are in the GTF
        gtf_tx_set = set(df_txlengths.transcript_id)
        tx_set = set(transcript_ids)

        # Find matching and missing IDs
        matching_ids = tx_set & gtf_tx_set
        missing_ids = tx_set - gtf_tx_set

        # If none match, raise an error
        if len(matching_ids) == 0:
            logging.error("ERROR: None of the provided transcript IDs are found in the GTF.")
            raise ValueError(
                "ERROR: None of the provided transcript IDs are found in the GTF. "
                "Please ensure the transcript list matches the GTF, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )

        # If some are missing
        if len(missing_ids) > 0:
            logging.error(f"{len(missing_ids)} transcript IDs not found in the GTF: {sorted(missing_ids)[:10]}{' ...' if len(missing_ids) > 10 else ''}")
            raise ValueError(
                "Some user-provided transcript IDs are missing from the GTF. "
                "Please ensure the transcript list matches the GTF, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )
        # Filter based on user provided transcript ID
        df_txlengths_filtered = df_txlengths.loc[df_txlengths.transcript_id.isin(transcript_ids)]


    else:
        # Select longest transcript per gene: hierarchy: CDS length > exon length; remaining ties are resolved alphabetically on transcript ID
        # Sort
        df_txlengths = df_txlengths.sort_values(by=['cds_length', 'exon_length', 'unspliced_length', 'transcript_id'], ascending=[False, False, False, True])
        # Drop duplicates on gene_id, keep first row
        df_txlengths_filtered = df_txlengths.drop_duplicates(subset='gene_id', keep='first')
        transcript_ids = df_txlengths_filtered.transcript_id.unique().tolist()
        transcript_ids_sorted = df_gtf.loc[df_gtf['Feature'] == 'transcript', 'transcript_id'].drop_duplicates().loc[lambda x: x.isin(df_txlengths_filtered.transcript_id)].tolist()



    # Check if all filtered genes have exactly one matching transcript ID
    # Count number of transcripts per gene
    transcripts_per_gene = df_txlengths_filtered.groupby("gene_id").transcript_id.nunique()
    # Raise error if any gene has not exactly one transcript and print problematic genes
    multi_assigned_genes = transcripts_per_gene[transcripts_per_gene > 1].index.tolist()
    if not (transcripts_per_gene == 1).all():
        logging.error("Some genes do not have exactly one representative transcript. Offending genes (more than 1 transcript assigned): ", multi_assigned_genes)
        raise ValueError(
            "ERROR: Some genes do not have exactly one representative transcript. "
            "Please make sure you provide a single transcript ID per gene ID, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )

    set_filtgenes = set(df_txlengths_filtered.gene_id)
    logging.info(f"Remaining genes after filtering: {len(set_filtgenes)}")

    # Check that all genes from GTF have been assigned at least 1 transcript
    missing_genes = input_genes - set_filtgenes
    if len(missing_genes) > 0:
        missing_genes_list = sorted(missing_genes)
        top = missing_genes_list[:10]
        more = f" and {len(missing_genes_list) - 10} more" if len(missing_genes_list) > 10 else ""

        logging.error(f"Not all genes in the GTF have a representative transcript. Genes without trancripts, after filtering: {top}{more}")
        raise ValueError(
            "ERROR: Not all genes in the GTF have a representative transcript after filtering. "
            "Please make sure you provide a single transcript ID for every gene ID in the GTF, "
            "or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )

    # Save transcript IDs to file
    logging.info(f"Saving representative transcript IDs to {output}.txt")
    with open(output + ".txt", "w") as f:
        f.write("\n".join(map(str, transcript_ids_sorted)) + "\n")

    # Make a transcript fai file
    # Set transcript_id as index
    df_txlengths_filtered.set_index('transcript_id', inplace=True)
    # Sort in order of transcript_ids
    df_txlengths_filtered = df_txlengths_filtered.loc[transcript_ids_sorted]
    # Reset index
    df_txlengths_filtered.reset_index(inplace=True)

    # Save transcript id and length (sum of all exons) to .fai file, without header
    logging.info(f"Saving representative transcript fai to {output}.fai")
    df_txlengths_filtered[['transcript_id', 'exon_length']].to_csv(f"{output}.fai", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

    # Create the transcript GTF file
    df_txlengths_filtered['src'] = 'src'
    df_txlengths_filtered['gene'] = 'gene'
    df_txlengths_filtered['start'] = 1
    df_txlengths_filtered['score'] = '.'
    df_txlengths_filtered['strand'] = '+'
    df_txlengths_filtered['frame'] = '.'
    df_txlengths_filtered['attributes'] = df_txlengths_filtered.apply(lambda row: f'gene_id:{row["gene_id"]}; transcript_id:{row["transcript_id"]}', axis=1)

    logging.info(f"Saving representative transcript GTF to {output}.gtf")
    # Order into gtf format
    df_txlengths_filtered[['transcript_id', 'src', 'gene', 'start', 'exon_length', 'score', 'strand', 'frame', 'attributes']].to_csv(f"{output}.gtf", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
    logging.info("Completed.")
    return


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--transcript", default="!{transcript}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.transcript, args.output)