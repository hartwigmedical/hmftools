package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionGeneFilters
{
    private final IsofoxConfig mConfig;
    private final List<ChrBaseRegion> mRestrictedGeneRegions;
    private final List<ChrBaseRegion> mExcludedGeneRegions;

    public FusionGeneFilters(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mExcludedGeneRegions = Lists.newArrayList();
        mRestrictedGeneRegions = Lists.newArrayList();

        buildGeneRegions(geneTransCache);
    }

    public boolean skipRead(final String otherChromosome, int otherPosition)
    {
        if(!HumanChromosome.contains(otherChromosome))
            return true;

        if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(otherChromosome))
            return true;
        else if(!mConfig.SpecificRegions.isEmpty() && !mConfig.SpecificRegions.stream().anyMatch(x -> x.containsPosition(otherChromosome, otherPosition)))
            return true;

        if(!mRestrictedGeneRegions.isEmpty() && !mRestrictedGeneRegions.stream().anyMatch(x -> x.containsPosition(otherChromosome, otherPosition)))
        {
            return true;
        }

        if(mExcludedGeneRegions.stream().anyMatch(x -> x.containsPosition(otherChromosome, otherPosition)))
        {
            return true;
        }

        return false;
    }

    private void buildGeneRegions(final EnsemblDataCache geneTransCache)
    {
        mConfig.EnrichedGeneIds.stream()
                .map(x -> geneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mExcludedGeneRegions.add(new ChrBaseRegion(
                        x.Chromosome, x.GeneStart - ENRICHED_GENE_BUFFER, x.GeneEnd + ENRICHED_GENE_BUFFER)));

        mConfig.ExcludedGeneIds.stream()
                .map(x -> geneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mExcludedGeneRegions.add(new ChrBaseRegion(
                        x.Chromosome, x.GeneStart, x.GeneEnd)));

        mConfig.RestrictedGeneIds.stream()
                .map(x -> geneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mRestrictedGeneRegions.add(new ChrBaseRegion(
                        x.Chromosome, x.GeneStart - 1000, x.GeneEnd + 1000)));
    }
}
