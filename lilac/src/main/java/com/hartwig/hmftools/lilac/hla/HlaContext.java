package com.hartwig.hmftools.lilac.hla;

import java.util.List;

public class HlaContext
{
    public final String Gene;
    public final List<Integer> AminoAcidBoundaries;
    public final com.hartwig.hmftools.lilac.fragment.ExpectedAlleles ExpectedAlleles;

    public HlaContext(final String gene, final List<Integer> aminoAcidBoundaries, final com.hartwig.hmftools.lilac.fragment.ExpectedAlleles expectedAlleles)
    {
        Gene = gene;
        AminoAcidBoundaries = aminoAcidBoundaries;
        ExpectedAlleles = expectedAlleles;
    }

}
