package com.hartwig.hmftools.redux.splice;

// NH (number-of-hits) -> MAPQ mapping. Encapsulated in a single class so the formula is easy to swap
// when comparing against aligners that use a different scheme. The current implementation matches STAR's
// default ladder (STAR ReadAlign.cpp): NH=1 -> 255, NH=2 -> 3, NH=3 -> 2, NH in [4,9] -> 1, NH>=10 -> 0.
//
// The ladder is meaningful only when NH reflects a normalized alignment set — duplicate-locus alts
// collapsed and discriminator-dropped losers excluded. Without that normalization, tx-contig
// representation artifacts inflate NH and the ladder over-penalizes legitimately-unique reads. See
// SpliceLiftBackConfig.STAR_MAPQ_LADDER_DESC for the recommended pairing.
public final class NhMapqLadder
{
    private NhMapqLadder() {}

    public static int starLadder(final int nh)
    {
        if(nh <= 1)
            return 255;
        if(nh == 2)
            return 3;
        if(nh == 3)
            return 2;
        if(nh < 10)
            return 1;
        return 0;
    }
}
