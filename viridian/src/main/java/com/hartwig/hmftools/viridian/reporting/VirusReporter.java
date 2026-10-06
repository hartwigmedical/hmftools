package com.hartwig.hmftools.viridian.reporting;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_COPIES_PER_CELL_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_INTEGRATIONS_MIN;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.OncologyGroupInfo;
import com.hartwig.hmftools.viridian.reference.VirusReference;

import org.jetbrains.annotations.Nullable;

// Decides if each virus is reported for downstream analysis and Orange report.
public class VirusReporter
{
    // One report per oncology group that is present or has an integration.
    public static List<VirusReport> report(
            List<DetectedVirus> detectedViruses, Map<OncologyGroup, Integer> integrationCounts,
            VirusReference reference, @Nullable Double expectedViralDepthPerCopy)
    {
        Map<OncologyGroup, DetectedVirus> presentByGroup = detectedViruses.stream()
                .collect(toMap(DetectedVirus::oncologyGroup, identity()));

        Set<OncologyGroup> groups = new HashSet<>(presentByGroup.keySet());
        groups.addAll(integrationCounts.keySet());

        return groups.stream()
                .map(group -> reportGroup(
                        group, presentByGroup.get(group), integrationCounts.getOrDefault(group, 0),
                        reference.oncologyGroupInfo(group), expectedViralDepthPerCopy))
                .toList();
    }

    static VirusReport reportGroup(
            OncologyGroup group, @Nullable DetectedVirus detected, int integrations,
            OncologyGroupInfo info, @Nullable Double expectedViralDepthPerCopy)
    {
        ContigStats representativeStats = detected != null ? detected.representativeContigStats() : null;

        // Null when it couldn't be calculated from the Purple data, in which case the virus can only be reported on integrations.
        Double copiesPerTumorCell = representativeStats != null && expectedViralDepthPerCopy != null
                ? representativeStats.depth().mean() / expectedViralDepthPerCopy
                : null;

        VirusReportStatus reason = decide(info.isReportable(), integrations, copiesPerTumorCell);

        return new VirusReport(
                group, info.reportingType(), info.driverLikelihood(), integrations,
                detected != null, representativeStats, copiesPerTumorCell, reason);
    }

    static VirusReportStatus decide(boolean reportable, int integrations, @Nullable Double copiesPerTumorCell)
    {
        if(!reportable)
        {
            return VirusReportStatus.NOT_REPORTABLE;
        }
        else if(integrations >= REPORTED_INTEGRATIONS_MIN)
        {
            return VirusReportStatus.REPORTED_ON_INTEGRATION;
        }
        else if(copiesPerTumorCell == null)
        {
            return VirusReportStatus.COPY_NUMBER_UNEVALUABLE;
        }
        else if(copiesPerTumorCell >= REPORTED_COPIES_PER_CELL_MIN)
        {
            return VirusReportStatus.REPORTED_ON_COPY_NUMBER;
        }
        else
        {
            return VirusReportStatus.COPY_NUMBER_TOO_LOW;
        }
    }
}
