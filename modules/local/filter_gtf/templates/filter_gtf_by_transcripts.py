#!/usr/bin/env python3

"""
Filters genonme annotation (GTF) by a list of selected transcript IDs (.txt).

Outputs:
- A filtered genome annotation GTF containing only entries corresponding to representative transcripts
"""

import platform
import argparse
from sys import exit
import pyranges as pr
import pandas as pd
import warnings
import logging


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

    if transcript != "":
        logging.info(f"Using the transcript IDs in: {transcript}")
        with open(transcript, "r") as file:
            transcript_ids = [line.strip() for line in file]
            # Remove potential empty strings caused by trailing newline
            transcript_ids = [x for x in transcript_ids if x]

        # Check if the transcript file was empty or only contained blank lines
        if len(transcript_ids) == 0:
            logging.error("Transcript file is empty.")
            raise ValueError(
                "ERROR: The transcript list provided is empty. "
                "Please provide a non-empty file with transcript IDs."
            )

        # Load the genome annotation
        df_gtf = parse_gtf(gtf)

        # Check if all user-provided IDs are in the GTF
        gtf_tx_set = set(df_gtf.transcript_id)
        tx_set = set(transcript_ids)

        # Find matching and missing IDs
        matching_ids = tx_set & gtf_tx_set
        missing_ids = tx_set - gtf_tx_set

        # If none match, raise an error
        if len(matching_ids) == 0:
            logging.error("None of the provided transcript IDs are found in the GTF.")
            raise ValueError(
                "ERROR: None of the provided transcript IDs are found in the GTF. "
                "Please ensure the transcript list matches the GTF, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )

        # If some are missing
        if len(missing_ids) > 0:
            logging.error(f"{len(missing_ids)} transcript IDs not found in the GTF: {sorted(missing_ids)[:10]}{' ...' if len(missing_ids) > 10 else ''}")
            raise ValueError(
                "ERROR: Some user-provided transcript IDs are missing from the GTF. "
                "Please ensure the transcript list matches the GTF, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )

        # Filter the genome annotation
        filt_annot = df_gtf.loc[(df_gtf['transcript_id'].isin(transcript_ids)) | (df_gtf['Feature'] == 'gene')].copy()


        # Count number of transcripts per gene in the filtered annotation
        transcripts_per_gene = filt_annot.groupby("gene_id").transcript_id.nunique()
        # Raise warning if any gene has not exactly one transcript and print problematic genes
        multi_assigned_genes = transcripts_per_gene[transcripts_per_gene > 1].index.tolist()
        if not (transcripts_per_gene == 1).all():
            logging.error("Some genes do not have exactly one representative transcript. Offending genes (more than 1 transcript assigned): ", multi_assigned_genes)
            raise ValueError(
                "ERROR: Some genes do not have exactly one representative transcript. "
                "Please make sure you provide a single transcript ID per gene ID, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )

        set_filtgenes = set(filt_annot.gene_id)
        logging.info(f"Remaining genes after filtering: {len(set_filtgenes)}")

        # Check that all genes from GTF have been assigned at least 1 transcript
        missing_genes = set(df_gtf.gene_id) - set_filtgenes
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


        # Convert the filtered annotation to pyranges
        pr_gtf = pr.PyRanges(filt_annot[[c for c in filt_annot.columns if c!='length']])

        # Save the filtered annotation
        logging.info(f"Saving genome GTF filtered by representative transcripts to {output}_filtered.gtf")
        pr_gtf.to_gtf(f"{output}_filtered.gtf")
        logging.info("Completed.")
    else:
        logging.error("No transcripts file found.")
        raise ValueError(
                "ERROR: No transcripts file found. "
                "Please provide a transcript list with transcript IDs that match the GTF, or omit --representative_transcript to let the pipeline automatically select representative transcripts."
            )

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


