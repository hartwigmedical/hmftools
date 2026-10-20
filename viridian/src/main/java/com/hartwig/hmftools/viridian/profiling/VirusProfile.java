package com.hartwig.hmftools.viridian.profiling;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeGenome;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeStatus;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// What the reads show about one oncology group present in the sample.
public record VirusProfile(
        GenomeStats representativeStats,
        RepresentativeStatus representativeStatus
)
{
    private static final Logger LOGGER = LogManager.getLogger(VirusProfile.class);

    public OncologyGroup oncologyGroup()
    {
        return representativeStats.genome().oncologyGroup();
    }

    public static List<VirusProfile> from(
            List<RepresentativeGenome> representatives, Map<VirusGenome, GenomeStats> representativeStats)
    {
        List<VirusProfile> profiles = new ArrayList<>();
        for(RepresentativeGenome representative : representatives)
        {
            GenomeStats stats = representativeStats.get(representative.genome());
            // Rarely, read assignment gives every read of a representative to other viruses' representatives, so it has
            // no measurement and its oncology group is dropped.
            if(stats == null)
            {
                LOGGER.warn(
                        "Representative genome {} assigned no reads; dropping oncology group {}",
                        representative.genome().name(), representative.oncologyGroup());
            }
            else
            {
                profiles.add(new VirusProfile(stats, representative.status()));
            }
        }
        return profiles;
    }
}
