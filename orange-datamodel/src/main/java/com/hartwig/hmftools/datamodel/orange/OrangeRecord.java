package com.hartwig.hmftools.datamodel.orange;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRecord
{
    @NotNull
    String sampleId();

    @NotNull
    LocalDate samplingDate();

    @NotNull
    ExperimentType experimentType();

    @NotNull
    Set<OrangeDoidNode> configuredPrimaryTumor();

    @NotNull
    OrangeRefGenomeVersion refGenomeVersion();

    @Nullable
    String platinumVersion();

    @Nullable
    OrangeSample refSample();

    @NotNull
    OrangeSample tumorSample();

    @Nullable
    Map<String, Double> germlineMVLHPerGene();

    @NotNull
    PurpleRecord purple();

    @NotNull
    LinxRecord linx();

    @NotNull
    List<WildTypeGene> wildTypeGenes();

    @Nullable
    IsofoxRecord isofox();

    @NotNull
    LilacRecord lilac();

    @NotNull
    ImmuneEscapeRecord immuneEscape();

    @Nullable
    VirusInterpreterData virusInterpreter();

    @Nullable
    ChordRecord chord();

    @Nullable
    CuppaData cuppa();

    @Nullable
    Set<PeachGenotype> peach();

    @Nullable
    List<SignatureAllocation> sigAllocations();

    @NotNull
    Map<PercentileType, Evaluation> cohortEvaluations();

    @NotNull
    OrangePlots plots();

    default boolean tumorOnlyMode()
    {
        return refSample() == null;
    }
}
