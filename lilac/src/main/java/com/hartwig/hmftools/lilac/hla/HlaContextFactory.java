package com.hartwig.hmftools.lilac.hla;

import java.util.List;

import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;

public class HlaContextFactory
{
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment;
    private final List<Integer> ABoundaries;
    private final List<Integer> BBoundaries;
    private final List<Integer> CBoundaries;

    public HlaContextFactory(
            final List<Integer> aBoundaries, final List<Integer> bBoundaries, final List<Integer> cBoundaries)
    {
        NucleotideGeneEnrichment = new NucleotideGeneEnrichment(aBoundaries, bBoundaries, cBoundaries);
        ABoundaries = aBoundaries;
        BBoundaries = bBoundaries;
        CBoundaries = cBoundaries;
    }

    public final HlaContext hlaA()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(NucleotideGeneEnrichment.getAFilterB(), NucleotideGeneEnrichment.getAFilterC());
        return new HlaContext("A", ABoundaries, expectedAlleles);
    }

    public final HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(NucleotideGeneEnrichment.getBFilterA(), NucleotideGeneEnrichment.getBFilterC());
        return new HlaContext("B", BBoundaries, expectedAlleles);
    }

    public final HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(NucleotideGeneEnrichment.getCFilterA(), NucleotideGeneEnrichment.getCFilterB());
        return new HlaContext("C", CBoundaries, expectedAlleles);
    }
}
