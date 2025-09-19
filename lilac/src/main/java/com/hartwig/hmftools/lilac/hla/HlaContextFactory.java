package com.hartwig.hmftools.lilac.hla;

import java.util.List;
import java.util.Map;
import java.util.NavigableMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;

public class HlaContextFactory
{
    private final LilacConfig mConfig;
    private final NucleotideGeneEnrichment NucleotideGeneEnrichment__;
    private final Map<HlaGene_, List<Integer>> GeneBoundaries_;

    public HlaContextFactory(final LilacConfig config, final Map<HlaGene_, List<Integer>> geneBoundaries_)
    {
        mConfig = config;
        NucleotideGeneEnrichment__ = NucleotideGeneEnrichment.create(geneBoundaries_);
        GeneBoundaries_ = geneBoundaries_;
    }

    public List<HlaContext> contexts()
    {
        List<HlaContext> output = Lists.newArrayList();
        int geneCount = mConfig.Genes.genes_().size();
        for(HlaGene_ gene : mConfig.Genes.genes_())
        {
            NavigableMap<Integer, Integer> filters = NucleotideGeneEnrichment__ == null ? Maps.newTreeMap() : NucleotideGeneEnrichment__.getFilters(gene);
            ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(geneCount, filters);
            HlaContext context = new HlaContext(gene, GeneBoundaries_.get(gene), expectedAlleles);
            output.add(context);
        }

        return output;
    }
}
