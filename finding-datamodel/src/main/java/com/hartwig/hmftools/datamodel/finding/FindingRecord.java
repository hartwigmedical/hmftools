package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record FindingRecord(
        @NotNull MetaProperties metaProperties,
        @NotNull PurityPloidyFit purityPloidyFit,
        @NotNull DriverFindingList<SmallVariant> smallVariants,
        @NotNull DriverFindingList<GainDeletion> gainDeletions,
        @NotNull DriverFindingList<Fusion> fusions,
        @NotNull DriverFindingList<Disruption> disruptions,
        @NotNull DriverFindingList<Virus> viruses,
        @NotNull FindingList<HlaAllele> hla,
        @NotNull FindingList<PharmocoGenotype> pharmocoGenotypes,
        @NotNull FindingItem<MicrosatelliteStability> microsatelliteStability,
        @NotNull FindingItem<TumorMutationStatus> tumorMutationStatus,
        @NotNull FindingItem<PredictedTumorOrigin> predictedTumorOrigin,
        @NotNull FindingItem<HomologousRecombination> homologousRecombination)
{}
