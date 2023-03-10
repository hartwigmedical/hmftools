package com.hartwig.hmftools.orange.algo;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.peach.PeachRecord;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
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
    public abstract Set<DoidNode> configuredPrimaryTumor();

    @NotNull
    public abstract RefGenomeVersion refGenomeVersion();

    @Nullable
    public abstract String platinumVersion();

    @Nullable
    public abstract OrangeSample refSample();

    @NotNull
    public abstract OrangeSample tumorSample();

    @Nullable
    public abstract Map<String, Double> germlineMVLHPerGene();

    @NotNull
    public abstract PurpleInterpretedData purple();

    @NotNull
    public abstract LinxRecord linx();

    @NotNull
    public abstract List<WildTypeGene> wildTypeGenes();

    @Nullable
    public abstract IsofoxInterpretedData isofox();

    @NotNull
    public abstract LilacSummaryData lilac();

    @Nullable
    public abstract VirusInterpreterData virusInterpreter();

    @Nullable
    public abstract ChordData chord();

    @Nullable
    public abstract CuppaData cuppa();

    @Nullable
    public abstract PeachRecord peach();

    @Nullable
    public abstract List<SignatureAllocation> sigAllocations();

    @NotNull
    public abstract Map<PercentileType, Evaluation> cohortEvaluations();

    @NotNull
    public abstract OrangePlots plots();

}
