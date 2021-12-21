package com.hartwig.hmftools.sage.coverage;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class Coverage
{
    private final Map<String,List<GeneCoverage>> mGeneCoverage;

    public Coverage(final Set<String> samples, final Collection<NamedBed> panel)
    {
        mGeneCoverage = Maps.newHashMap();
        final Set<String> genes = panel.stream().map(NamedBed::name).collect(Collectors.toSet());

        samples.forEach(sample -> mGeneCoverage.put(sample, createGeneCoverage(genes, panel)));
    }

    public Set<String> samples()
    {
        return mGeneCoverage.keySet();
    }

    public List<GeneDepth> depth(final String sample)
    {
        return coverage(sample).stream()
                .map(GeneCoverage::geneDepth)
                .sorted(Comparator.comparing(GeneDepth::gene))
                .collect(Collectors.toList());
    }

    public List<GeneCoverage> coverage(final String sample)
    {
        return mGeneCoverage.getOrDefault(sample, Collections.emptyList());
    }

    public List<GeneCoverage> coverage(final String sample, final String chromosome)
    {
        return coverage(sample).stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList());
    }

    private List<GeneCoverage> createGeneCoverage(final Set<String> genes, final Collection<NamedBed> panel)
    {
        final List<GeneCoverage> result = Lists.newArrayList();
        for(String gene : genes)
        {
            List<NamedBed> exons = panel.stream().filter(x -> x.name().equals(gene)).collect(Collectors.toList());
            if(!exons.isEmpty())
            {
                result.add(new GeneCoverage(exons));
            }
        }

        return result;
    }
}
