package com.hartwig.hmftools.lilac.hla;

import java.util.Set;

import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles;

public class HlaContext
{
    public final String Gene;
    public final Set<Integer> AminoAcidBoundaries;
    public final com.hartwig.hmftools.lilac.nuc.ExpectedAlleles ExpectedAlleles;

    public HlaContext(final String gene, final Set<Integer> aminoAcidBoundaries, final com.hartwig.hmftools.lilac.nuc.ExpectedAlleles expectedAlleles)
    {
        Gene = gene;
        AminoAcidBoundaries = aminoAcidBoundaries;
        ExpectedAlleles = expectedAlleles;
    }

}
