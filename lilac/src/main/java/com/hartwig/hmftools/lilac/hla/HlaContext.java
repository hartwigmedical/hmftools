package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;

import java.util.List;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;

public class HlaContext
{
    public final String Gene;
    public final List<Integer> AminoAcidBoundaries;
    public final ExpectedAlleles ExpectedAlleles;

    public HlaContext(final String gene, final List<Integer> aminoAcidBoundaries, final ExpectedAlleles expectedAlleles)
    {
        Gene = gene;
        AminoAcidBoundaries = aminoAcidBoundaries;
        ExpectedAlleles = expectedAlleles;
    }

    public String geneName() { return HLA_PREFIX + Gene; }

}
