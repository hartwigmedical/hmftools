package com.hartwig.hmftools.redux.splice;

public final class SpliceCommon
{
    public static final String CONTIG_NAME_DELIM = "_";
    public static final String CONTIG_NAME_PREFIX = "ens";

    public static final String TRANSCRIPT_CONTIGS_FILE_ID = ".fasta";
    public static final String CONTIG_MAPPINGS_FILE_ID = ".rna_contigs_mappings.tsv";

    public static final String ALT_CONTIG_SUFFIX = "_tx";

    public static String altContigName(final String chromosome)
    {
        return chromosome + ALT_CONTIG_SUFFIX;
    }
}
