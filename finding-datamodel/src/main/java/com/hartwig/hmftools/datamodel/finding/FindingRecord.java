package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

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
    @NotNull
    DriverFindings<SmallVariant> smallVariants();

    @NotNull
    DriverFindings<GainDeletion> gainDeletions();

    @NotNull
    DriverFindings<Fusion> fusions();

    @NotNull
    DriverFindings<Disruption> disruptions();

    @NotNull
    DriverFindings<Virus> viruses();

    @NotNull
    Findings<HlaAllele> hla();

    @NotNull
    Findings<PharmocoGenotype> pharmocoGenotypes();

    @NotNull
    List<ChromosomeArmCopyNumber> chromosomeArmCopyNumbers();

    @NotNull CharacteristicsFinding<MicrosatelliteStability> microsatelliteStability();

    @NotNull CharacteristicsFinding<TumorMutationStatus> tumorMutationStatus();

    @NotNull CharacteristicsFinding<PredictedTumorOrigin> predictedTumorOrigin();

    @NotNull CharacteristicsFinding<HomologousRecombination> homologousRecombination();

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
