package com.hartwig.hmftools.viridian.reporting;

import static java.util.Objects.requireNonNull;
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

import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.reference.VirusGroup;
import com.hartwig.hmftools.viridian.reference.VirusGroupInfo;

import org.jetbrains.annotations.Nullable;

// Decides if each virus is reported for downstream analysis and Orange report.
public class VirusReporter
{
    // One report per virus group that is profiled or has an integration.
    public static List<VirusReport> report(
            List<VirusProfile> profiles, Map<VirusGroup, Integer> integrationCounts,
            Function<VirusGroup, VirusGroupInfo> virusGroupInfo, ClonalityInputs clonalityInputs)
    {
        Map<VirusGroup, VirusProfile> profileByGroup = profiles.stream()
                .collect(toMap(VirusProfile::virusGroup, identity()));

        Set<VirusGroup> groupsToReport = Stream.concat(profileByGroup.keySet().stream(), integrationCounts.keySet().stream())
                .collect(Collectors.toSet());

        return groupsToReport.stream()
                .map(group -> reportGroup(
                        group, profileByGroup.get(group), integrationCounts.getOrDefault(group, 0),
                        virusGroupInfo.apply(group), clonalityInputs))
                .toList();
    }

    static VirusReport reportGroup(
            VirusGroup group, @Nullable VirusProfile profile, int integrations, VirusGroupInfo info,
            ClonalityInputs clonalityInputs)
    {
        Double copiesPerTumorCell = null;
        ClonalityStatus clonalityStatus;
        if(profile == null)
        {
            clonalityStatus = ClonalityStatus.NOT_PROFILED;
        }
        else if(clonalityInputs.purpleFailure() == null)
        {
            copiesPerTumorCell = profile.representativeStats().depth().mean() / requireNonNull(clonalityInputs.expectedViralDepthPerCopy());
            clonalityStatus = copiesPerTumorCell >= REPORTED_COPIES_PER_CELL_MIN ? ClonalityStatus.CLONAL : ClonalityStatus.NOT_CLONAL;
        }
        else
        {
            clonalityStatus = switch(clonalityInputs.purpleFailure())
            {
                case FAIL_NO_TUMOR -> ClonalityStatus.PURPLE_FAIL_NO_TUMOR;
                case FAIL_CONTAMINATION -> ClonalityStatus.PURPLE_FAIL_CONTAMINATION;
                default -> throw new IllegalStateException(
                        "Not a clonality-blocking Purple failure: " + clonalityInputs.purpleFailure());
            };
        }

        boolean hasMinIntegrations = integrations >= REPORTED_INTEGRATIONS_MIN;
        boolean isReported = info.isReportable() && (hasMinIntegrations || clonalityStatus == ClonalityStatus.CLONAL);

        return new VirusReport(
                group, info.reportingType(), integrations, hasMinIntegrations, profile,
                copiesPerTumorCell, clonalityStatus, isReported);
    }
}
