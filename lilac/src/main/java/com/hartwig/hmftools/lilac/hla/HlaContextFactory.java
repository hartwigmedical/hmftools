package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_1;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_2;
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
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment;
    private final Map<HlaGene, List<Integer>> GeneBoundaries_;

    public HlaContextFactory(final Map<HlaGene, List<Integer>> geneBoundaries_)
    {
        NucleotideGeneEnrichment = new NucleotideGeneEnrichment(geneBoundaries_);
        GeneBoundaries_ = geneBoundaries_;
    }

    private HlaContext hlaA()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getAFilterB_(), NucleotideGeneEnrichment.getAFilterC_());

        return new HlaContext(HLA_A, GeneBoundaries_.get(HLA_A), expectedAlleles);
    }

    private HlaContext hlaB()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getBFilterA_(), NucleotideGeneEnrichment.getBFilterC_());

        return new HlaContext(HLA_B, GeneBoundaries_.get(HLA_B), expectedAlleles);
    }

    private HlaContext hlaC()
    {
        ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(
                NucleotideGeneEnrichment.getCFilterA_(), NucleotideGeneEnrichment.getCFilterB_());

        return new HlaContext(HLA_C, GeneBoundaries_.get(HLA_C), expectedAlleles);
    }

    public List<HlaContext> contexts(final LilacConfig config)
    {
        List<HlaContext> output = Lists.newArrayList();
        if(config.ClassType == null || config.ClassType == CLASS_1)
        {
            output.add(hlaA());
            output.add(hlaB());
            output.add(hlaC());
        }

        if(config.ClassType == null || config.ClassType == CLASS_2)
        {
            for(HlaGene gene : HlaGene.values())
            {
                if(gene.isPseudo() || gene.mhcClass() != CLASS_2)
                    continue;

                output.add(new HlaContext(gene, GeneBoundaries_.get(gene), new ExpectedAlleles(null)));
            }
        }

        return output;
    }
}
