package com.hartwig.hmftools.sage;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_38;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bamops.BamSampler;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SageCommon
{
    public static final String APP_NAME = "Sage";

    public static final String SAMPLE_DELIM = ",";

    public static final Logger SG_LOGGER = LogManager.getLogger(SageCommon.class);

    public static void setReadLength(
            final SageConfig config, final Map<Chromosome,List<BaseRegion>> panelRegions, final String bamFile)
    {
        if(config.getReadLength() > 0) // skip if set in config
            return;

        BamSampler bamSampler = new BamSampler(config.RefGenomeFile);

        List<ChrBaseRegion> sampleRegions = Lists.newArrayList();

        if(!config.SpecificChrRegions.Regions.isEmpty())
        {
            sampleRegions.addAll(config.SpecificChrRegions.Regions);
        }
        else if(!panelRegions.isEmpty())
        {
            RefGenomeCoordinates refGenomeCoordinates = config.RefGenVersion.is37() ? COORDS_37 : COORDS_38;

            for(Map.Entry<Chromosome, List<BaseRegion>> entry : panelRegions.entrySet())
            {
                String chromosome = config.RefGenVersion.versionedChromosome(entry.getKey().toString());
                int chromosomeLength = refGenomeCoordinates.length(chromosome);

                for(BaseRegion region : entry.getValue())
                {
                    // if the region is a single-base hotspot, build some width around it
                    int regionStart = region.start();
                    int regionEnd = region.end();

                    if(region.length() == 1)
                    {
                        regionStart = max(regionStart - 500, 0);
                        regionEnd = min(regionEnd + 500, chromosomeLength);
                    }

                    sampleRegions.add(new ChrBaseRegion(chromosome, regionStart, regionEnd));
                }
            }
        }
        else
        {
            sampleRegions.add(bamSampler.defaultRegion());
        }

        for(ChrBaseRegion sampleRegion : sampleRegions)
        {
            if(bamSampler.calcBamCharacteristics(bamFile, sampleRegion) && bamSampler.maxReadLength() > 0)
            {
                config.setReadLength(bamSampler.maxReadLength());
                return;
            }
        }

        SG_LOGGER.debug("BAM read-length sampling using default read length({})", DEFAULT_READ_LENGTH);
        config.setReadLength(DEFAULT_READ_LENGTH);
    }
}
