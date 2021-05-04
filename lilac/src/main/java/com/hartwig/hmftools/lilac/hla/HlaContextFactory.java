package com.hartwig.hmftools.lilac.hla;

import java.util.Set;

import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.nuc.NucleotideGeneEnrichment;

public class HlaContextFactory
{
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment;
    private final Set<Integer> ABoundaries;
    private final Set<Integer> BBoundaries;
    private final Set<Integer> CBoundaries;

    public HlaContextFactory(
            final Set<Integer> aBoundaries, final Set<Integer> bBoundaries, final Set<Integer> cBoundaries)
    {
        NucleotideGeneEnrichment = new NucleotideGeneEnrichment(aBoundaries, bBoundaries, cBoundaries);
        ABoundaries = aBoundaries;
        BBoundaries = bBoundaries;
        CBoundaries = cBoundaries;
    }

    public final HlaContext hlaA()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(this.NucleotideGeneEnrichment.getAFilterB(), this.NucleotideGeneEnrichment.getAFilterC());
        return new HlaContext("A", this.ABoundaries, expectedAlleles);
    }

    public final HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(this.NucleotideGeneEnrichment.getBFilterA(), this.NucleotideGeneEnrichment.getBFilterC());
        return new HlaContext("B", this.BBoundaries, expectedAlleles);
    }

    public final HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(this.NucleotideGeneEnrichment.getCFilterA(), this.NucleotideGeneEnrichment.getCFilterB());
        return new HlaContext("C", this.CBoundaries, expectedAlleles);
    }
}
