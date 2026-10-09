package com.hartwig.hmftools.viridian.detection;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupResolution;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeContigCandidate;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.jetbrains.annotations.Nullable;

// Holds the final detection status info for each detected oncology group.
public record DetectedVirus(
        OncologyGroup oncologyGroup,
        OncologyGroupResolution resolution,
        OncologyGroupOutcome outcome,
        // Reads aligning to any genome in the group.
        int groupReadCount,
        // Genomes in this group which had any alignments.
        int alignedContigCount,
        // Genomes that passed the support filter.
        int candidateCount,
        // Genomes abundant enough to contest the lead. A measure of how ambiguous the group is.
        int comparableCandidateCount,
        // Null where no representative was resolved.
        @Nullable ContigStats representativeContigStats
)
{
    public DetectedVirus
    {
        if(resolution == OncologyGroupResolution.NO_CANDIDATES)
        {
            throw new IllegalArgumentException("Oncology group is not detected: " + oncologyGroup);
        }
        if((resolution == OncologyGroupResolution.RESOLVED) == (representativeContigStats == null))
        {
            throw new IllegalArgumentException("Representative stats do not match resolution: " + resolution);
        }
    }

    public static List<DetectedVirus> from(
            List<OncologyGroupRepresentativeSelection> selections, Map<ViralContig, ContigStats> representativeStats,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        return selections.stream()
                .filter(selection -> selection.resolution() != OncologyGroupResolution.NO_CANDIDATES)
                .map(selection -> from(selection, representativeStats, groupReadCounts))
                .toList();
    }

    private static DetectedVirus from(
            OncologyGroupRepresentativeSelection selection, Map<ViralContig, ContigStats> representativeStats,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        OncologyGroup oncologyGroup = selection.oncologyGroup();
        int comparableCandidates = (int) selection.candidates().stream()
                .filter(RepresentativeContigCandidate::comparable).count();

        ViralContig representative = selection.representative();
        ContigStats contigStats = representative != null ? representativeStats.get(representative) : null;
        if(representative != null && contigStats == null)
        {
            throw new IllegalStateException("Representative was not measured: " + representative.name());
        }

        return new DetectedVirus(
                oncologyGroup, selection.resolution(), selection.outcome(),
                groupReadCounts.getOrDefault(oncologyGroup, 0),
                selection.candidates().size() + selection.rejected().size(),
                selection.candidates().size(), comparableCandidates, contigStats);
    }
}
