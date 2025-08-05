package com.hartwig.hmftools.lilac.hla;

import java.util.List;

import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;

import org.apache.commons.lang3.NotImplementedException;

public class HlaContext
{
    public final HlaGene Gene;
    public final List<Integer> AminoAcidBoundaries;
    public final ExpectedAlleles ExpectedAlleles;

    public HlaContext(final HlaGene gene, final List<Integer> aminoAcidBoundaries, final ExpectedAlleles expectedAlleles)
    {
        // TODO:
        if(aminoAcidBoundaries == null)
        {
            throw new RuntimeException("aminoAcidBoundaries == null");
        }

        // TODO:
        if(expectedAlleles == null)
        {
            throw new RuntimeException("expectedAlleles == null");
        }

        Gene = gene;
        AminoAcidBoundaries = aminoAcidBoundaries;
        ExpectedAlleles = expectedAlleles;
    }

    public String geneName() { return Gene.toString(); }
}
