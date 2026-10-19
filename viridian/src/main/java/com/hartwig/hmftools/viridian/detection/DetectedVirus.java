package com.hartwig.hmftools.viridian.detection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeGenome;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeStatus;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Holds the final status info for each detected oncology group.
public record DetectedVirus(
        GenomeStats representativeStats,
        RepresentativeStatus representativeStatus
)
{
    private static final Logger LOGGER = LogManager.getLogger(DetectedVirus.class);

    public OncologyGroup oncologyGroup()
    {
        return representativeStats.genome().oncologyGroup();
    }

    public static List<DetectedVirus> from(
            List<RepresentativeGenome> representatives, Map<VirusGenome, GenomeStats> representativeStats)
    {
        List<DetectedVirus> detectedViruses = new ArrayList<>();
        for(RepresentativeGenome representative : representatives)
        {
            GenomeStats stats = representativeStats.get(representative.genome());
            // Rarely, read assignment gives every read of a representative to other viruses' representatives, so it has
            // no measurement and its oncology group is not detected.
            if(stats == null)
            {
                LOGGER.warn(
                        "Representative genome {} assigned no reads; dropping oncology group {}",
                        representative.genome().name(), representative.oncologyGroup());
            }
            else
            {
                detectedViruses.add(new DetectedVirus(stats, representative.status()));
            }
        }
        return detectedViruses;
    }
}
