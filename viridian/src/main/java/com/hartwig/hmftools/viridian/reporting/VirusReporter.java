package com.hartwig.hmftools.viridian.reporting;

import static java.util.Objects.requireNonNull;
import static java.util.function.Function.identity;
import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;
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
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;
import com.hartwig.hmftools.viridian.reference.VirusGroupInfo;

import org.jetbrains.annotations.Nullable;

// Decides if each virus is reported for downstream analysis and Orange report.
public class VirusReporter
{
    // One report per virus group that is profiled or has an integration.
    public static List<VirusReport> report(
            List<VirusProfile> profiles, List<VirusGenome> integratedGenomes,
            Function<VirusGroup, VirusGroupInfo> virusGroupInfo, ClonalityInputs clonalityInputs)
    {
        Map<VirusGroup, VirusProfile> profileByGroup = profiles.stream()
                .collect(toMap(VirusProfile::virusGroup, identity()));
        Map<VirusGroup, List<VirusGenome>> integratedGenomesByGroup = integratedGenomes.stream()
                .collect(groupingBy(VirusGenome::virusGroup));

        Set<VirusGroup> groupsToReport = Stream.concat(profileByGroup.keySet().stream(), integratedGenomesByGroup.keySet().stream())
                .collect(Collectors.toSet());

        return groupsToReport.stream()
                .map(group -> reportGroup(
                        group, profileByGroup.get(group), integratedGenomesByGroup.getOrDefault(group, List.of()),
                        virusGroupInfo.apply(group), clonalityInputs))
                .toList();
    }

    static VirusReport reportGroup(
            VirusGroup group, @Nullable VirusProfile profile, List<VirusGenome> integratedGenomes, VirusGroupInfo info,
            ClonalityInputs clonalityInputs)
    {
        Double copiesPerTumorCell = null;
        ClonalityStatus clonalityStatus;
        if(profile == null)
        {
            copiesPerTumorCell = clonalityInputs.purpleFailure() == null ? 0.0 : null;
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

        int integrations = integratedGenomes.size();
        boolean isReported = info.isReportable()
                && (integrations >= REPORTED_INTEGRATIONS_MIN || clonalityStatus == ClonalityStatus.CLONAL);

        List<String> integratedVirusNames = integratedGenomes.stream().map(VirusGenome::virusName).toList();
        String virusName = profile != null
                ? profile.representativeStats().genome().virusName()
                // The most common virus name from the integrations.
                : integratedVirusNames.stream()
                        .collect(groupingBy(identity(), counting()))
                        .entrySet().stream()
                        .min(Map.Entry.<String, Long>comparingByValue().reversed().thenComparing(Map.Entry.comparingByKey()))
                        .orElseThrow()
                        .getKey();
        ReportQcStatus qcStatus = integratedVirusNames.stream().allMatch(virusName::equals)
                ? ReportQcStatus.NO_ABNORMALITIES
                : ReportQcStatus.AMBIGUOUS_VIRUS;

        return new VirusReport(
                group, virusName, qcStatus, info.reportingType(), integrations, profile, copiesPerTumorCell,
                clonalityStatus, isReported);
    }
}
