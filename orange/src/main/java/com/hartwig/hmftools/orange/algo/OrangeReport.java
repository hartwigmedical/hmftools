package com.hartwig.hmftools.orange.algo;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.*;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.Set;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangeReport {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract LocalDate experimentDate();

    @NotNull
    public abstract ExperimentType experimentType();

    @NotNull
    public abstract Set<OrangeDoidNode> configuredPrimaryTumor();

    @NotNull
    public abstract OrangeRefGenomeVersion refGenomeVersion();

    @Nullable
    public abstract String platinumVersion();

    @Nullable
    public abstract OrangeSample refSample();

    @NotNull
    public abstract OrangeSample tumorSample();

    @Nullable
    public abstract Map<String, Double> germlineMVLHPerGene();

    @NotNull
    public abstract PurpleRecord purple();

    @NotNull
    public abstract LinxRecord linx();

    @NotNull
    public abstract List<WildTypeGene> wildTypeGenes();

    @Nullable
    public abstract IsofoxInterpretedData isofox();

    @NotNull
    public abstract LilacRecord lilac();

    @Nullable
    public abstract VirusInterpreterData virusInterpreter();

    @Nullable
    public abstract ChordRecord chord();

    @Nullable
    public abstract CuppaData cuppa();

    @Nullable
    public abstract Set<PeachGenotype> peach();

    @Nullable
    public abstract List<SignatureAllocation> sigAllocations();

    @NotNull
    public abstract Map<PercentileType, Evaluation> cohortEvaluations();

    @NotNull
    public abstract OrangePlots plots();

}
