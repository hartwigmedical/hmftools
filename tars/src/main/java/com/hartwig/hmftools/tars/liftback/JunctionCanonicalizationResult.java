package com.hartwig.hmftools.tars.liftback;

// Result of JunctionCanonicalizer.tryCanonicalize. changed=false means no junction was shifted
// (newCigar is unset). The alignment start never moves: only interior M|N|M boundaries are
// re-partitioned, so the leftmost mapped base is unaffected.
public record JunctionCanonicalizationResult(boolean changed, String newCigar)
{
    public static JunctionCanonicalizationResult unchanged()
    {
        return new JunctionCanonicalizationResult(false, null);
    }
}
