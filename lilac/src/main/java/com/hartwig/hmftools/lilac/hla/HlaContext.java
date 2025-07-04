package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.GeneCache.longGeneName;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;

public class HlaContext
{
    public final String Gene;
    public final Set<Integer> AminoAcidBoundaries;
    public final ExpectedAlleles ExpectedAlleles;

    public HlaContext(final String gene, final Iterable<Integer> aminoAcidBoundaries, final ExpectedAlleles expectedAlleles)
    {
        Gene = gene;
        AminoAcidBoundaries = Sets.newTreeSet(aminoAcidBoundaries);
        ExpectedAlleles = expectedAlleles;
    }

    public String geneName() { return longGeneName(Gene); }

}
