package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Stream;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.IFindingList;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingRecord(
        @NotNull String version,
        @NotNull MetaProperties metaProperties,
        @NotNull Qc qc,
        @NotNull PurityPloidyFit purityPloidyFit,
        @NotNull DriverFindingList<SmallVariant> somaticSmallVariants,
        @NotNull DriverFindingList<SmallVariant> germlineSmallVariants,
        @NotNull DriverFindingList<GainDeletion> somaticGainDeletions,
        @NotNull DriverFindingList<GainDeletion> germlineGainDeletions,
        @NotNull DriverFindingList<Disruption> somaticDisruptions,
        @NotNull DriverFindingList<Disruption> germlineDisruptions,
        @NotNull DriverFindingList<Fusion> fusions,
        @NotNull DriverFindingList<Virus> viruses,
        @NotNull FindingList<ChromosomeArmCopyNumber> chromosomeArmCopyNumbers,
        @NotNull FindingList<HlaAllele> hlaAlleles,
        @NotNull FindingList<PharmacoGenotype> pharmacoGenotypes,
        @NotNull FindingItem<PredictedTumorOrigin> predictedTumorOrigin,
        @NotNull FindingItem<MicrosatelliteStability> microsatelliteStability,
        @NotNull FindingItem<TumorMutationalLoad> tumorMutationalLoad,
        @NotNull FindingItem<TumorMutationalBurden> tumorMutationalBurden,
        @NotNull FindingItem<HomologousRecombination> homologousRecombination)
{
    public static final String VERSION = FindingRecord.class.getPackage().getImplementationVersion();

    public boolean hasLowPurity()
    {
        List<IFindingList<?>> lists =
                List.of(somaticSmallVariants, somaticGainDeletions, somaticDisruptions, fusions, viruses, chromosomeArmCopyNumbers, hlaAlleles, pharmacoGenotypes);
        for(IFindingList<?> list : lists)
        {
            if(list.status().errors().contains(FindingStatus.Issue.LOW_PURITY))
            {
                return true;
            }
            if(list.status().warnings().contains(FindingStatus.Issue.LOW_PURITY))
            {
                return true;
            }
        }
        return false;
    }

    @NotNull
    public DriverFindingList<SmallVariant> allSmallVariants()
    {
        return mergeDriverFindings(somaticSmallVariants, germlineSmallVariants);
    }

    @NotNull
    public DriverFindingList<GainDeletion> allGainDeletions()
    {
        return mergeDriverFindings(somaticGainDeletions, germlineGainDeletions);
    }

    @NotNull
    public DriverFindingList<Disruption> allDisruptions()
    {
        return mergeDriverFindings(somaticDisruptions, germlineDisruptions);
    }

    @NotNull
    public List<SmallVariant> allReportableSmallVariants()
    {
        return allSmallVariants().reportedOnly().findings();
    }

    @NotNull
    public List<GainDeletion> allReportableGainDeletions()
    {
        return allGainDeletions().reportedOnly().findings();
    }

    @NotNull
    public List<Disruption> allReportableDisruptions()
    {
        return allDisruptions().reportedOnly().findings();
    }

    private static <T extends Driver> DriverFindingList<T> mergeDriverFindings(DriverFindingList<T> somatic, DriverFindingList<T> germline)
    {
        List<T> merged = mergeAndSort(somatic.findings(), germline.findings());
        return DriverFindingListBuilder.<T>builder()
                .status(mergeFindingStatus(somatic.status(), germline.status()))
                .findings(merged).build();
    }

    @NotNull
    private static <T extends Driver> List<T> mergeAndSort(List<T> list1, List<T> list2)
    {
        return Stream.concat(list1.stream(), list2.stream())
                .sorted(T.COMPARATOR)
                .toList();
    }

    static FindingStatus mergeFindingStatus(FindingStatus somaticStatus, FindingStatus germlineStatus)
    {

        SortedSet<FindingStatus.Issue> errors = new TreeSet<>();
        SortedSet<FindingStatus.Issue> warnings = new TreeSet<>();
        warnings.addAll(somaticStatus.warnings());
        warnings.addAll(germlineStatus.warnings());
        FindingStatus.Status status = switch(somaticStatus.status())
        {
            case OK -> switch(germlineStatus.status())
            {
                case OK:
                    yield FindingStatus.Status.OK;
                case NOT_RELIABLE, NOT_AVAILABLE:
                    warnings.addAll(germlineStatus.errors());
                    yield FindingStatus.Status.OK;
            };
            case NOT_RELIABLE -> switch(germlineStatus.status())
            {
                case OK:
                    warnings.addAll(somaticStatus.errors());
                    yield FindingStatus.Status.OK;
                case NOT_RELIABLE, NOT_AVAILABLE:
                    yield FindingStatus.Status.NOT_RELIABLE;
            };
            case NOT_AVAILABLE -> switch(germlineStatus.status())
            {
                case OK:
                    warnings.addAll(somaticStatus.errors());
                    yield FindingStatus.Status.OK;
                case NOT_RELIABLE:
                    yield FindingStatus.Status.NOT_RELIABLE;
                case NOT_AVAILABLE:
                    yield FindingStatus.Status.NOT_AVAILABLE;
            };
        };
        if(status != FindingStatus.Status.OK)
        {
            errors.addAll(somaticStatus.errors());
            errors.addAll(germlineStatus.errors());
        }
        return FindingStatusBuilder.builder()
                .status(status)
                .errors(errors)
                .warnings(warnings)
                .build();
    }
}
