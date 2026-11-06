package com.hartwig.hmftools.viridian.reporting;

import static java.util.Objects.requireNonNull;
import static java.util.function.Function.identity;
import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_COPIES_PER_CELL_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_INTEGRATIONS_MIN;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.reference.Virus;
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
                        profileByGroup.get(group), integratedGenomesByGroup.getOrDefault(group, List.of()),
                        virusGroupInfo.apply(group), clonalityInputs))
                .toList();
    }

    static VirusReport reportGroup(
            @Nullable VirusProfile profile, List<VirusGenome> integratedGenomes, VirusGroupInfo info,
            ClonalityInputs clonalityInputs)
    {
        Clonality clonality = assessClonality(profile, clonalityInputs);

        int integrations = integratedGenomes.size();
        boolean isReported = info.isReportable()
                && (integrations >= REPORTED_INTEGRATIONS_MIN || clonality.status() == ClonalityStatus.CLONAL);

        List<Virus> integratedViruses = integratedGenomes.stream().map(VirusGenome::virus).toList();
        // Within the virus group, want to report the specific virus too.
        // Usually read profiling and integrations detect the same virus, so show that.
        // In rare cases where they don't agree, prefer the profiling and show a warning.
        Virus virus = profile != null ? profile.representativeStats().genome().virus() : mostCommonVirus(integratedViruses);
        ReportQcStatus qcStatus = integratedViruses.stream().allMatch(virus::equals)
                ? ReportQcStatus.NO_ABNORMALITIES
                : ReportQcStatus.AMBIGUOUS_VIRUS;

        return new VirusReport(
                virus, qcStatus, info.reportingType(), integrations, profile, clonality.copiesPerTumorCell(),
                clonality.status(), isReported);
    }

    private record Clonality(
            @Nullable Double copiesPerTumorCell,
            ClonalityStatus status
    )
    {
    }

    private static Clonality assessClonality(@Nullable VirusProfile profile, ClonalityInputs inputs)
    {
        if(profile == null)
        {
            return new Clonality(inputs.purpleFailure() == null ? 0.0 : null, ClonalityStatus.NOT_PROFILED);
        }
        else if(inputs.purpleFailure() == null)
        {
            double copiesPerTumorCell = profile.representativeStats().depth().mean() / requireNonNull(inputs.expectedViralDepthPerCopy());
            ClonalityStatus status =
                    copiesPerTumorCell >= REPORTED_COPIES_PER_CELL_MIN ? ClonalityStatus.CLONAL : ClonalityStatus.NOT_CLONAL;
            return new Clonality(copiesPerTumorCell, status);
        }
        else
        {
            ClonalityStatus status = switch(inputs.purpleFailure())
            {
                case FAIL_NO_TUMOR -> ClonalityStatus.PURPLE_FAIL_NO_TUMOR;
                case FAIL_CONTAMINATION -> ClonalityStatus.PURPLE_FAIL_CONTAMINATION;
                default -> throw new IllegalStateException(
                        "Not a clonality-blocking Purple failure: " + inputs.purpleFailure());
            };
            return new Clonality(null, status);
        }
    }

    private static Virus mostCommonVirus(List<Virus> viruses)
    {
        Map<Virus, Long> counts = viruses.stream().collect(groupingBy(identity(), counting()));
        Comparator<Map.Entry<Virus, Long>> mostCommonFirst =
                Map.Entry.<Virus, Long>comparingByValue().reversed().thenComparing(Map.Entry.comparingByKey());
        return counts.entrySet().stream().min(mostCommonFirst).orElseThrow().getKey();
    }
}
