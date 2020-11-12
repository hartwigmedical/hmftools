package com.hartwig.hmftools.sage.coverage;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class Coverage {

    private final Map<String, List<GeneCoverage>> geneCoverage;

    public Coverage(final Collection<String> referenceSamples, final Collection<String> tumorSamples, final Collection<NamedBed> panel) {
        this.geneCoverage = Maps.newHashMap();
        final Set<String> genes = panel.stream().map(NamedBed::name).collect(Collectors.toSet());

        tumorSamples.forEach(sample -> geneCoverage.put(sample, createGeneCoverage(genes, panel)));
        referenceSamples.forEach(sample -> geneCoverage.put(sample, createGeneCoverage(genes, panel)));
    }

    @NotNull
    public Set<String> samples() {
        return geneCoverage.keySet();
    }

    @NotNull
    public List<GeneDepth> depth(@NotNull final String sample) {
        return coverage(sample).stream().map(x -> x.geneDepth(30)).collect(Collectors.toList());
    }

    @NotNull
    public List<GeneCoverage> coverage(@NotNull final String sample) {
        return geneCoverage.getOrDefault(sample, Collections.emptyList());
    }

    @NotNull
    public List<GeneCoverage> coverage(@NotNull final String sample, @NotNull final String chromosome) {
        return coverage(sample).stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList());
    }

    @NotNull
    private List<GeneCoverage> createGeneCoverage(Set<String> genes, Collection<NamedBed> panel) {
        final List<GeneCoverage> result = Lists.newArrayList();
        for (String gene : genes) {
            List<NamedBed> exons = panel.stream().filter(x -> x.name().equals(gene)).collect(Collectors.toList());
            if (!exons.isEmpty()) {
                result.add(new GeneCoverage(exons));
            }
        }

        return result;
    }

}
