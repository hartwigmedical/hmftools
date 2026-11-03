package com.hartwig.hmftools.viridian.profiling;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeSelection;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// What the reads show about one oncology group present in the sample.
public record VirusProfile(
        GenomeStats representativeStats,
        List<VirusGenome> secondaryGenomes
)
{
    private static final Logger LOGGER = LogManager.getLogger(VirusProfile.class);

    public OncologyGroup oncologyGroup()
    {
        return representativeStats.genome().oncologyGroup();
    }

    public static List<VirusProfile> from(
            List<RepresentativeSelection> selections, Map<VirusGenome, GenomeStats> representativeStats)
    {
        List<VirusProfile> profiles = new ArrayList<>();
        for(RepresentativeSelection selection : selections)
        {
            VirusGenome representative = selection.representativeGenome();
            GenomeStats stats = representativeStats.get(representative);
            // Rarely, read assignment gives every read of a representative to other viruses' representatives, so it has
            // no measurement and its oncology group is dropped.
            if(stats == null)
            {
                LOGGER.warn(
                        "Representative genome {} assigned no reads; dropping oncology group {}",
                        representative.contigName(), representative.oncologyGroup());
            }
            else
            {
                profiles.add(new VirusProfile(stats, selection.secondaryGenomes()));
            }
        }
        return profiles;
    }
}
