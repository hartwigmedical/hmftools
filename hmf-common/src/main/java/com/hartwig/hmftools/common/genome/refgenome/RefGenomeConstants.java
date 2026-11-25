package com.hartwig.hmftools.common.genome.refgenome;

public final class RefGenomeConstants
{
    public static final double GENOME_LENGTH_V37 = 3.088e10;
    public static final double GENOME_LENGTH_V38 = 3.096e10;

    public static final long BASES_PER_GENOME = 2_859_000_000L;
    public static final double MB_PER_GENOME = BASES_PER_GENOME * 1e-6;

    public static final long CODING_BASES_PER_GENOME = 31_880_000; // calculated from GRCh38 canonical transcripts (overlaps ignored)
}
