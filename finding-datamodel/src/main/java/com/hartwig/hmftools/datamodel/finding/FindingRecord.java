package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FindingRecord
{
    @NotNull DriverFindingList<SmallVariant> smallVariants();

    @NotNull DriverFindingList<GainDeletion> gainDeletions();

    @NotNull DriverFindingList<Fusion> fusions();

    @NotNull DriverFindingList<Disruption> disruptions();

    @NotNull DriverFindingList<Virus> viruses();

    @NotNull FindingList<HlaAllele> hla();

    @NotNull FindingList<PharmocoGenotype> pharmocoGenotypes();

    @NotNull FindingItem<MicrosatelliteStability> microsatelliteStability();

    @NotNull FindingItem<TumorMutationStatus> tumorMutationStatus();

    @NotNull FindingItem<PredictedTumorOrigin> predictedTumorOrigin();

    @NotNull FindingItem<HomologousRecombination> homologousRecombination();

    @NotNull
    OrangeRefGenomeVersion refGenomeVersion();

    @NotNull
    ExperimentType experimentType();

    @Nullable
    String pipelineVersion();

    @NotNull
    String version();

    @NotNull
    PurpleFit purpleFit();
}
