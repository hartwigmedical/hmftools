package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.region.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class TargetRegions
{
    private final Map<String,List<BaseRegion>> mTargetRegions;

    public TargetRegions(final String targetRegionsBedFile, final RefGenomeVersion refGenVersion)
    {
        mTargetRegions = Maps.newHashMap();

        if(targetRegionsBedFile != null)
        {
            loadTargetRegionsBed(targetRegionsBedFile, refGenVersion);
        }
    }

    public boolean hasTargetRegions() { return !mTargetRegions.isEmpty(); }

    public boolean inTargetRegions(final String chromsome, int position)
    {
        final List<BaseRegion> chrRegions = mTargetRegions.get(chromsome);

        if(chrRegions == null)
            return false;

        return chrRegions.stream().anyMatch(x -> positionWithin(position, x.start(), x.end()));
    }

    private void loadTargetRegionsBed(final String filename, final RefGenomeVersion refGenVersion)
    {
        if(filename == null)
            return;

        Map<Chromosome,List<BaseRegion>> chrRegionsMap = loadBedFileChrMap(filename);

        for(Map.Entry<Chromosome,List<BaseRegion>> entry : chrRegionsMap.entrySet())
        {
            String chromosome = refGenVersion.versionedChromosome(entry.getKey().toString());
            mTargetRegions.put(chromosome, entry.getValue());

        }

        SV_LOGGER.info("loaded {} target regions from file({})",
                mTargetRegions.values().stream().mapToInt(x -> x.size()).sum(), filename);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(TARGET_REGIONS_BED, false, TARGET_REGIONS_BED_DESC);
    }
}
