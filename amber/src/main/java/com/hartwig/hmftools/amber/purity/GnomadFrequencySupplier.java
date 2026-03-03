package com.hartwig.hmftools.amber.purity;

import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public interface GnomadFrequencySupplier
{
    double getFrequency(String chromosome, int position);
}

class DefaultGnomadFrequencySupplier implements GnomadFrequencySupplier
{
    private final Map<String, Map<Integer, Double>> mGnomadData = new HashMap<>();

    public DefaultGnomadFrequencySupplier(ListMultimap<Chromosome, AmberSite> mChromosomeSites, RefGenomeVersion genomeVersion)
    {
        for(Chromosome chromosome : mChromosomeSites.keySet())
        {
            Map<Integer, Double> chromosomeData = new HashMap<>();
            for(AmberSite site : mChromosomeSites.get(chromosome))
            {
                chromosomeData.put(site.position(), site.VariantAlleleFrequency);
            }
            mGnomadData.put(genomeVersion.versionedChromosome(chromosome), chromosomeData);
        }
    }

    @Override
    public double getFrequency(String chromosome, int position)
    {
        if(mGnomadData.containsKey(chromosome))
        {
            Map<Integer, Double> chromosomeData = mGnomadData.get(chromosome);
            return chromosomeData.getOrDefault(position, 0.0);
        }
        throw new IllegalArgumentException("No Gnomad data for: " + chromosome + ":" + position);
    }
}
