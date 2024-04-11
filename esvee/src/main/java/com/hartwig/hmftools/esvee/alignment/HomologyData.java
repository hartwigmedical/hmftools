package com.hartwig.hmftools.esvee.alignment;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class HomologyData
{
    public final String Homology;
    public final int ExactStart;
    public final int ExactEnd;
    public final int InexactStart;
    public final int InexactEnd;

    public HomologyData(final String homology, final int exactStart, final int exactEnd, final int inexactStart, final int inexactEnd)
    {
        Homology = homology;
        ExactStart = exactStart;
        ExactEnd = exactEnd;
        InexactStart = inexactStart;
        InexactEnd = inexactEnd;
    }

    public static HomologyData determineHomology(
            final AlignData alignmentStart, final AlignData alignmentEnd, final RefGenomeInterface refGenome)
    {

        return null;
    }
}
