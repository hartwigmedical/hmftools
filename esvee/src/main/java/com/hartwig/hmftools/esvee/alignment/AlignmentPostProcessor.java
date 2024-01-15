package com.hartwig.hmftools.esvee.alignment;

import com.hartwig.hmftools.esvee.old.AlignedAssembly;
import com.hartwig.hmftools.esvee.old.Alignment;

public enum AlignmentPostProcessor
{
    ;

    /** Flips an assembly if it contains more inverted bases than non-inverted, or if it contains an equal number but
     * starts inverted */
    public static AlignedAssembly flipIfRequired(final AlignedAssembly assembly)
    {
        int mapped = 0;
        int mappedInverted = 0;
        for(Alignment alignment : assembly.getAlignmentBlocks())
        {
            if(!alignment.isMapped())
                continue;
            mapped++;
            if(alignment.Inverted)
                mappedInverted++;
        }

        final int mappedNonInverted = mapped - mappedInverted;
        if(mappedInverted > mappedNonInverted)
            return assembly.flipStrand();
        return assembly;
    }
}
