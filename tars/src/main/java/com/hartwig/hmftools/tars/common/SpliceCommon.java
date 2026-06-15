package com.hartwig.hmftools.tars.common;

public final class SpliceCommon
{
    // Minimum M flank (bp) for an N to count as a trusted splice junction. A flank below this carries too
    // little evidence to assert a junction (~1/4^8 chance a coincidental anchor matches), so an N with a
    // sub-threshold flank is treated as fabricated/untrustworthy: the discriminator won't swap the primary
    // onto it, and the rescue-fold / terminal-collapse passes re-evaluate the terminal anchor against the
    // genome rather than trusting bwa's placement. One threshold, used everywhere a junction anchor is judged.
    public static final int MIN_JUNCTION_ANCHOR = 8;

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
