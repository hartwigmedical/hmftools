package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.MhcClass.CLASS_1;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;

public class HlaContextFactory
{
    private final LilacConfig mConfig;
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment;
    private final Map<HlaGene, List<Integer>> GeneBoundaries;

    public HlaContextFactory(final LilacConfig config, final Map<HlaGene, List<Integer>> geneBoundaries)
    {
        mConfig = config;
        NucleotideGeneEnrichment = config.Genes.coversMhcClass1() ? new NucleotideGeneEnrichment(geneBoundaries) : null;
        GeneBoundaries = geneBoundaries;
    }

    private HlaContext hlaA()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getAFilterB(), NucleotideGeneEnrichment.getAFilterC());

        return new HlaContext(HLA_A, GeneBoundaries.get(HLA_A), expectedAlleles);
    }

    private HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getBFilterA(), NucleotideGeneEnrichment.getBFilterC());

        return new HlaContext(HLA_B, GeneBoundaries.get(HLA_B), expectedAlleles);
    }

    private HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getCFilterA(), NucleotideGeneEnrichment.getCFilterB());

        return new HlaContext(HLA_C, GeneBoundaries.get(HLA_C), expectedAlleles);
    }

    public List<HlaContext> contexts()
    {
        List<HlaContext> output = Lists.newArrayList();
        if(mConfig.Genes.coversMhcClass1())
        {
            output.add(hlaA());
            output.add(hlaB());
            output.add(hlaC());
        }

        for(HlaGene gene : mConfig.Genes.genes())
        {
            if(gene.isPseudo() || gene.mhcClass() == CLASS_1)
                continue;

            output.add(new HlaContext(gene, GeneBoundaries.get(gene), new ExpectedAlleles(null)));
        }

        return output;
    }
}
