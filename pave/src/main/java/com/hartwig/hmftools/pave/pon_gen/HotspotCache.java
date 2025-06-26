package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.pon_gen.PonConfig.GERMLINE_HOTSPOT;
import static com.hartwig.hmftools.pave.pon_gen.PonConfig.SOMATIC_HOTSPOT;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;

public class HotspotCache
{
    private final ListMultimap<Chromosome, VariantHotspot> mSomaticHotspotMap;
    private final ListMultimap<Chromosome, VariantHotspot> mGermlineHotspotMap;

    public HotspotCache(final ConfigBuilder configBuilder)
    {
        mSomaticHotspotMap = ArrayListMultimap.create();
        mGermlineHotspotMap = ArrayListMultimap.create();

        try
        {
            if(configBuilder.hasValue(SOMATIC_HOTSPOT))
            {
                String hotspotFile = configBuilder.getValue(SOMATIC_HOTSPOT);
                mSomaticHotspotMap.putAll(VariantHotspotFile.readFromVCF(hotspotFile));
            }
        }
        catch(Exception e)
        {
            PV_LOGGER.error("error loading hotspot file: {}", e.toString());
            System.exit(1);
        }

        try
        {
            if(configBuilder.hasValue(GERMLINE_HOTSPOT))
            {
                String hotspotFile = configBuilder.getValue(GERMLINE_HOTSPOT);
                mGermlineHotspotMap.putAll(VariantHotspotFile.readFromVCF(hotspotFile));
            }
        }
        catch(Exception e)
        {
            PV_LOGGER.error("error loading hotspot file: {}", e.toString());
            System.exit(1);
        }

        if(!mSomaticHotspotMap.isEmpty() || !mGermlineHotspotMap.isEmpty())
        {
            PV_LOGGER.info("loaded hotspots somatic({}) germline({})", mSomaticHotspotMap.size(), mGermlineHotspotMap.size());
        }
    }

    public List<VariantHotspot> getChromosomeSomaticHotspots(final String chromosome)
    {
        return mSomaticHotspotMap.get(HumanChromosome.fromString(chromosome));
    }

    public List<VariantHotspot> getChromosomeGermlineHotspots(final String chromosome)
    {
        return mGermlineHotspotMap.get(HumanChromosome.fromString(chromosome));
    }
}
