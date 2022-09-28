package com.hartwig.hmftools.sage.coverage;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

public class Coverage
{
    private final Map<String,List<GeneCoverage>> mGeneCoverage;

    private static final String GENE_EXON_DELIM = "_";

    public Coverage(final List<String> samples, final Collection<NamedBed> panel)
    {
        mGeneCoverage = Maps.newHashMap();

        if(panel.isEmpty())
            return;

        final Set<String> genes = panel.stream().map(x -> extractGeneName(x.name())).collect(Collectors.toSet());

        List<GeneCoverage> geneCoverages = Lists.newArrayList();

        for(String gene : genes)
        {
            List<NamedBed> exons = panel.stream().filter(x -> extractGeneName(x.name()).equals(gene)).collect(Collectors.toList());

            if(!exons.isEmpty())
            {
                List<ExonCoverage> exonCoverages = Lists.newArrayList();

                for(int i = 0; i < exons.size(); ++i)
                {
                    NamedBed exon = exons.get(i);
                    int exonRank = extractExonRank(exon.name(), i + 1);
                    exonCoverages.add(new ExonCoverage(exon, exonRank));
                }

                GeneCoverage geneCoverage = new GeneCoverage(gene, exonCoverages);
                geneCoverages.add(geneCoverage);
            }
        }

        mGeneCoverage.put(samples.get(0), geneCoverages);

        for(int i = 1; i < samples.size(); ++i)
        {
            // copy since depth is unique to the sample
            List<GeneCoverage> sampleGeneCoverages = geneCoverages.stream().map(x -> x.clone()).collect(Collectors.toList());
            mGeneCoverage.put(samples.get(i), sampleGeneCoverages);
        }
    }

    public Set<String> samples() { return mGeneCoverage.keySet(); }

    public void writeFiles(final String outputVcf)
    {
        if(mGeneCoverage.isEmpty())
            return;

        String parent = new File(outputVcf).getParent();

        for(Map.Entry<String,List<GeneCoverage>> entry : mGeneCoverage.entrySet())
        {
            String sampleId = entry.getKey();
            List<GeneCoverage> geneCoverages = entry.getValue();

            String geneFile = parent + File.separator + sampleId + ".sage.gene.coverage.tsv";
            String exonFile = parent + File.separator + sampleId + ".sage.exon.medians.tsv";

            GeneDepthFile.write(geneFile, depth(sampleId), GeneDepthBuilder.DEPTH_BUCKETS);
            ExonMedianDepth.write(exonFile, geneCoverages);
        }
    }

    private List<GeneDepth> depth(final String sample)
    {
        return coverage(sample).stream()
                .map(x -> GeneDepthBuilder.buildGeneDepth(x))
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

    private static String extractGeneName(final String name)
    {
        if(name.contains(GENE_EXON_DELIM))
            return name.split(GENE_EXON_DELIM, -1)[0];

        return name;
    }

    private static int extractExonRank(final String name, final int index)
    {
        if(name.contains(GENE_EXON_DELIM))
            return Integer.parseInt(name.split(GENE_EXON_DELIM, -1)[1]);

        return index;
    }
}
