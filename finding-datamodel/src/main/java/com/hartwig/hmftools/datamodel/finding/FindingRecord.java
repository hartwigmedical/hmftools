package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record FindingRecord(
        @NotNull ExperimentType experimentType,
        @Nullable String pipelineVersion,
        @NotNull String version,
        @NotNull OrangeRefGenomeVersion refGenomeVersion,
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
