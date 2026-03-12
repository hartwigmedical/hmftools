package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.stream.Stream;

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
        @NotNull FindingList<HlaAllele> hlaAlleles,
        @NotNull FindingList<PharmacoGenotype> pharmacoGenotypes,
        @NotNull FindingList<PredictedTumorOrigin> predictedTumorOrigins,
        @NotNull FindingItem<MicrosatelliteStability> microsatelliteStability,
        @NotNull FindingItem<TumorMutationalLoad> tumorMutationalLoad,
        @NotNull FindingItem<TumorMutationalBurden> tumorMutationalBurden,
        @NotNull FindingItem<HomologousRecombination> homologousRecombination,
        @NotNull VisualisationFiles visualisationFiles)
{
    @NotNull
    public List<SmallVariant> allSmallVariants()
    {
        return mergeDriverFindings(somaticSmallVariants, germlineSmallVariants);
    }

    @NotNull
    public List<GainDeletion> allGainDeletions()
    {
        return mergeDriverFindings(somaticGainDeletions, germlineGainDeletions);
    }

    @NotNull
    public List<Disruption> allDisruptions()
    {
        return mergeDriverFindings(somaticDisruptions, germlineDisruptions);
    }

    @NotNull
    public List<SmallVariant> allReportableSmallVariants()
    {
        return mergeReportedDriverFindings(somaticSmallVariants, germlineSmallVariants);
    }

    @NotNull
    public List<GainDeletion> allReportableGainDeletions()
    {
        return mergeReportedDriverFindings(somaticGainDeletions, germlineGainDeletions);
    }

    @NotNull
    public List<Disruption> allReportableDisruptions()
    {
        return mergeReportedDriverFindings(somaticDisruptions, germlineDisruptions);
    }

    @NotNull
    private static <T extends Driver> List<T> mergeDriverFindings(DriverFindingList<T> somatic, DriverFindingList<T> germline)
    {
        return mergeAndSort(somatic.findings(), germline.findings());
    }

    @NotNull
    private static <T extends Driver> List<T> mergeReportedDriverFindings(DriverFindingList<T> somatic, DriverFindingList<T> germline)
    {
        return mergeAndSort(somatic.reportedOnly().findings(), germline.reportedOnly().findings());
    }

    @NotNull
    private static <T extends Driver> List<T> mergeAndSort(List<T> list1, List<T> list2)
    {
        return Stream.concat(list1.stream(), list2.stream())
                .sorted(T.COMPARATOR)
                .toList();
    }
}
