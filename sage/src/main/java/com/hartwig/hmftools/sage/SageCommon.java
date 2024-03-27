package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSampler;

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

        ChrBaseRegion sampleRegion = null;

        if(!config.SpecificChrRegions.Regions.isEmpty())
        {
            sampleRegion = config.SpecificChrRegions.Regions.get(0);
        }
        else if(!panelRegions.isEmpty())
        {
            for(Map.Entry<Chromosome, List<BaseRegion>> entry : panelRegions.entrySet())
            {
                BaseRegion region = entry.getValue().get(0);

                sampleRegion = new ChrBaseRegion(
                        config.RefGenVersion.versionedChromosome(entry.getKey().toString()), region.start(), region.end());

                break;
            }
        }
        else
        {
            sampleRegion = bamSampler.defaultRegion();
        }

        if(bamSampler.calcBamCharacteristics(bamFile, sampleRegion) && bamSampler.maxReadLength() > 0)
        {
            config.setReadLength(bamSampler.maxReadLength());
        }
        else
        {
            SG_LOGGER.debug("BAM read-length sampling using default read length({})", DEFAULT_READ_LENGTH);
            config.setReadLength(DEFAULT_READ_LENGTH);
        }
    }
}
