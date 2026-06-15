package com.hartwig.hmftools.tars.liftback;

// Result of JunctionCanonicalizer.tryCanonicalize. Changed=false means no junction was shifted
// (NewCigar is unset). The alignment start never moves: only interior M|N|M boundaries are
// re-partitioned, so the leftmost mapped base is unaffected.
public class JunctionCanonicalizationResult
{
    public final boolean Changed;
    public final String NewCigar;
    public final int JunctionsShifted;

    public JunctionCanonicalizationResult(final boolean changed, final String newCigar, final int junctionsShifted)
    {
        Changed = changed;
        NewCigar = newCigar;
        JunctionsShifted = junctionsShifted;
    }

    public static JunctionCanonicalizationResult unchanged()
    {
        return new JunctionCanonicalizationResult(false, null, 0);
    }
}
