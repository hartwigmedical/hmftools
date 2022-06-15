package com.hartwig.hmftools.orange.algo;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.cuppa.CuppaData;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.lilac.LilacData;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangeReport {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract LocalDate reportDate();

    @NotNull
    public abstract Set<DoidNode> configuredPrimaryTumor();

    @NotNull
    public abstract RefGenomeVersion refGenomeVersion();

    @Nullable
    public abstract String platinumVersion();

    @NotNull
    public abstract OrangeSample refSample();

    @NotNull
    public abstract OrangeSample tumorSample();

    @NotNull
    public abstract Map<String, Double> germlineMVLHPerGene();

    @NotNull
    public abstract PurpleInterpretedData purple();

    @NotNull
    public abstract LinxInterpretedData linx();

    @Nullable
    public abstract IsofoxInterpretedData isofox();

    @NotNull
    public abstract LilacData lilac();

    @NotNull
    public abstract VirusInterpreterData virusInterpreter();

    @NotNull
    public abstract ChordAnalysis chord();

    @NotNull
    public abstract CuppaData cuppa();

    @NotNull
    public abstract List<PeachGenotype> peach();

    @NotNull
    public abstract List<ProtectEvidence> protect();

    @NotNull
    public abstract Map<PercentileType, Evaluation> cohortEvaluations();

    @NotNull
    public abstract OrangePlots plots();

}
