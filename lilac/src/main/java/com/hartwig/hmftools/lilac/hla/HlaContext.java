package com.hartwig.hmftools.lilac.hla;

import java.util.List;

import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;

public class HlaContext
{
    public final HlaGene_ Gene;
    public final List<Integer> AminoAcidBoundaries_;
    public final ExpectedAlleles ExpectedAlleles;

    public HlaContext(final HlaGene_ gene, final List<Integer> aminoAcidBoundaries_, final ExpectedAlleles expectedAlleles)
    {
        // TODO:
        if(aminoAcidBoundaries_ == null)
        {
            throw new RuntimeException("aminoAcidBoundaries == null");
        }

        // TODO:
        if(expectedAlleles == null)
        {
            throw new RuntimeException("expectedAlleles == null");
        }

        Gene = gene;
        AminoAcidBoundaries_ = aminoAcidBoundaries_;
        ExpectedAlleles = expectedAlleles;
    }

    public String geneName() { return Gene.toString(); }
}
