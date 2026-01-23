package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingRecord(
        @NotNull MetaProperties metaProperties,
        @NotNull PurityPloidyFit purityPloidyFit,
        @NotNull DriverFindingList<SmallVariant> somaticSmallVariants,
        @NotNull DriverFindingList<SmallVariant> germlineSmallVariants,
        @NotNull DriverFindingList<GainDeletion> somaticGainDeletions,
        @NotNull DriverFindingList<GainDeletion> germlineGainDeletions,
        @NotNull DriverFindingList<Disruption> somaticDisruptions,
        @NotNull DriverFindingList<Disruption> germlineDisruptions,
        @NotNull DriverFindingList<Fusion> fusions,
        @NotNull DriverFindingList<Virus> viruses,
        @NotNull FindingList<HlaAllele> hla,
        @NotNull FindingList<PharmocoGenotype> pharmocoGenotypes,
        @NotNull FindingItem<MicrosatelliteStability> microsatelliteStability,
        @NotNull FindingItem<TumorMutationStatus> tumorMutationStatus,
        @NotNull FindingItem<PredictedTumorOrigin> predictedTumorOrigin,
        @NotNull FindingItem<HomologousRecombination> homologousRecombination)
{
    @NotNull
    List<SmallVariant> allSmallVariants()
    {
        return mergeDriverFindings(somaticSmallVariants, germlineSmallVariants);
    }

    @NotNull
    List<GainDeletion> allGainDeletions()
    {
        return mergeDriverFindings(somaticGainDeletions, germlineGainDeletions);
    }

    @NotNull
    List<Disruption> allDisruptions()
    {
        return mergeDriverFindings(somaticDisruptions, germlineDisruptions);
    }

    @NotNull
    List<SmallVariant> allReportableSmallVariants()
    {
        return mergeReportedDriverFindings(somaticSmallVariants, germlineSmallVariants);
    }

    @NotNull
    List<GainDeletion> allReportableGainDeletions()
    {
        return mergeReportedDriverFindings(somaticGainDeletions, germlineGainDeletions);
    }

    @NotNull
    List<Disruption> allReportableDisruptions()
    {
        return mergeReportedDriverFindings(somaticDisruptions, germlineDisruptions);
    }

    @NotNull
    private static <T extends Driver> List<T> mergeDriverFindings(DriverFindingList<T> somatic, DriverFindingList<T> germline)
    {
        List<T> all = new ArrayList<>(somatic.findings());
        all.addAll(germline.findings());
        all.sort(T.COMPARATOR);
        return all;
    }

    @NotNull
    private static <T extends Driver> List<T> mergeReportedDriverFindings(DriverFindingList<T> somatic, DriverFindingList<T> germline)
    {
        List<T> all = new ArrayList<>(somatic.reportedOnly().findings());
        all.addAll(germline.reportedOnly().findings());
        all.sort(T.COMPARATOR);
        return all;
    }
}
