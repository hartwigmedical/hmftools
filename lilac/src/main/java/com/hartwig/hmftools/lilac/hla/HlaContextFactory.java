package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_1;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_2;
import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_C;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;

public class HlaContextFactory
{
    private final LilacConfig mConfig;
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment_;
    private final Map<HlaGene_, List<Integer>> GeneBoundaries_;

    public HlaContextFactory(final LilacConfig config, final Map<HlaGene_, List<Integer>> geneBoundaries_)
    {
        mConfig = config;
        NucleotideGeneEnrichment_ = config.Genes.coversMhcClass1() ? new NucleotideGeneEnrichment(geneBoundaries_) : null;
        GeneBoundaries_ = geneBoundaries_;
    }

    private HlaContext hlaA()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment_.getAFilterB_(), NucleotideGeneEnrichment_.getAFilterC_());

        return new HlaContext(HLA_A, GeneBoundaries_.get(HLA_A), expectedAlleles);
    }

    private HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment_.getBFilterA_(), NucleotideGeneEnrichment_.getBFilterC_());

        return new HlaContext(HLA_B, GeneBoundaries_.get(HLA_B), expectedAlleles);
    }

    private HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment_.getCFilterA_(), NucleotideGeneEnrichment_.getCFilterB_());

        return new HlaContext(HLA_C, GeneBoundaries_.get(HLA_C), expectedAlleles);
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

        for(HlaGene_ gene : mConfig.Genes.genes_())
        {
            if(gene.isPseudo() || gene.mhcClass() == CLASS_1)
                continue;

            output.add(new HlaContext(gene, GeneBoundaries_.get(gene), new ExpectedAlleles(null)));
        }

        return output;
    }
}
