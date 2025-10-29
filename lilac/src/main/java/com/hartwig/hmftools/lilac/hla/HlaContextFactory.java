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
    private final NucleotideGeneEnrichment mNucleotideGeneEnrichment;
    private final Map<HlaGene, List<Integer>> GeneBoundaries;

    public HlaContextFactory(final LilacConfig config, final Map<HlaGene, List<Integer>> geneBoundaries)
    {
        mConfig = config;
        mNucleotideGeneEnrichment = NucleotideGeneEnrichment.create(geneBoundaries);
        GeneBoundaries = geneBoundaries;
    }

    public List<HlaContext> contexts()
    {
        List<HlaContext> output = Lists.newArrayList();
        int geneCount = mConfig.Genes.genes().size();
        for(HlaGene gene : mConfig.Genes.genes())
        {
            NavigableMap<Integer, Integer> filters =
                    mNucleotideGeneEnrichment == null ? Maps.newTreeMap() : mNucleotideGeneEnrichment.getFilters(gene);
            ExpectedAlleles expectedAlleles = ExpectedAlleles.expectedAlleles(geneCount, filters);
            HlaContext context = new HlaContext(gene, GeneBoundaries.get(gene), expectedAlleles);
            output.add(context);
        }

        return output;
    }
}
