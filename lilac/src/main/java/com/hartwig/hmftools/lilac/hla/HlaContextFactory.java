package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;

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
        return new HlaContext(GENE_A, ABoundaries, expectedAlleles);
    }

    public final HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(NucleotideGeneEnrichment.getBFilterA(), NucleotideGeneEnrichment.getBFilterC());
        return new HlaContext(GENE_B, BBoundaries, expectedAlleles);
    }

    public final HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles =
                ExpectedAlleles.expectedAlleles(NucleotideGeneEnrichment.getCFilterA(), NucleotideGeneEnrichment.getCFilterB());
        return new HlaContext(GENE_C, CBoundaries, expectedAlleles);
    }
}
