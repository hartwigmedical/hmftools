package com.hartwig.hmftools.viridian.reporting;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_COPIES_PER_CELL_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_INTEGRATIONS_MIN;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.OncologyGroupInfo;

import org.jetbrains.annotations.Nullable;

// Decides if each virus is reported for downstream analysis and Orange report.
public class VirusReporter
{
    // One report per oncology group that is detected or has an integration.
    public static List<VirusReport> report(
            List<DetectedVirus> detectedViruses, Map<OncologyGroup, Integer> integrationCounts,
            Function<OncologyGroup, OncologyGroupInfo> oncologyGroupInfo, @Nullable Double expectedViralDepthPerCopy)
    {
        Map<OncologyGroup, DetectedVirus> detectedByGroup = detectedViruses.stream()
                .collect(toMap(DetectedVirus::oncologyGroup, identity()));

        Set<OncologyGroup> groupsToReport = Stream.concat(detectedByGroup.keySet().stream(), integrationCounts.keySet().stream())
                .collect(Collectors.toSet());

        return groupsToReport.stream()
                .map(group -> reportGroup(
                        group, detectedByGroup.get(group), integrationCounts.getOrDefault(group, 0),
                        oncologyGroupInfo.apply(group), expectedViralDepthPerCopy))
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

        VirusReportStatus status = decideStatus(info.isReportable(), integrations, copiesPerTumorCell);

        return new VirusReport(
                group, info.reportingType(), info.driverLikelihood(), integrations,
                detected != null, representativeStats, copiesPerTumorCell, status);
    }

    static VirusReportStatus decideStatus(boolean reportable, int integrations, @Nullable Double copiesPerTumorCell)
    {
        if(!reportable)
        {
            return VirusReportStatus.NOT_REPORTABLE;
        }

        boolean integrated = integrations >= REPORTED_INTEGRATIONS_MIN;
        boolean clonal = copiesPerTumorCell != null && copiesPerTumorCell >= REPORTED_COPIES_PER_CELL_MIN;

        if(integrated && clonal)
        {
            return VirusReportStatus.REPORTED_ON_INTEGRATION_AND_CLONALITY;
        }
        else if(integrated)
        {
            return VirusReportStatus.REPORTED_ON_INTEGRATION;
        }
        else if(clonal)
        {
            return VirusReportStatus.REPORTED_ON_CLONALITY;
        }
        else if(copiesPerTumorCell == null)
        {
            return VirusReportStatus.CLONALITY_UNEVALUABLE;
        }
        else
        {
            return VirusReportStatus.NOT_CLONAL;
        }
    }
}
