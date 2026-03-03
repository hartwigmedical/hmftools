package com.hartwig.hmftools.amber.purity;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
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
    private final Map<String, AmberSite[]> mGnomadData = new HashMap<>();

    public DefaultGnomadFrequencySupplier(ListMultimap<Chromosome, AmberSite> chromosomeSites, RefGenomeVersion genomeVersion)
    {
        for(Chromosome chromosome : chromosomeSites.keySet())
        {
            List<AmberSite> sites = chromosomeSites.get(chromosome);
            AmberSite[] sortedSites = sites.toArray(new AmberSite[0]);
            Arrays.sort(sortedSites, Comparator.comparingInt(AmberSite::position));
            mGnomadData.put(genomeVersion.versionedChromosome(chromosome), sortedSites);
        }
    }

    @Override
    public double getFrequency(String chromosome, int position)
    {
        AmberSite[] sites = mGnomadData.get(chromosome);
        if(sites != null)
        {
            int index = binarySearch(sites, position);
            if(index >= 0)
            {
                return sites[index].VariantAlleleFrequency;
            }
        }
        throw new IllegalArgumentException("No Gnomad data for: " + chromosome + ":" + position);
    }

    private static int binarySearch(AmberSite[] sites, int position)
    {
        int low = 0;
        int high = sites.length - 1;
        while(low <= high)
        {
            int mid = (low + high) >>> 1;
            int midPosition = sites[mid].Position;
            if(midPosition < position)
            {
                low = mid + 1;
            }
            else if(midPosition > position)
            {
                high = mid - 1;
            }
            else
            {
                return mid;
            }
        }
        return -1;
    }
}
