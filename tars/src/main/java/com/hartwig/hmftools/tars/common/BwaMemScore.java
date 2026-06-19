package com.hartwig.hmftools.tars.common;

// bwa-mem2 default affine alignment scoring, the single source for every liftback pass that re-scores a
// boundary: the BoundaryReclaim base walk (match/mismatch) and the LiftBackResolver CIGAR reconstruction
// (all four). Soft-clips are not penalised (clipped bases score 0).
public final class BwaMemScore
{
    public static final int MATCH = 1;
    public static final int MISMATCH = -4;
    public static final int GAP_OPEN = -6;
    public static final int GAP_EXTEND = -1;

    private BwaMemScore() { }
}
