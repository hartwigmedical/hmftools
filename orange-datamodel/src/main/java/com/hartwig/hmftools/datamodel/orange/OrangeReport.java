package com.hartwig.hmftools.datamodel.orange;

import com.hartwig.hmftools.datamodel.chord.ChordData;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.doid.DoidNode;
import com.hartwig.hmftools.datamodel.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.hla.LilacSummaryData;
import com.hartwig.hmftools.datamodel.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.datamodel.linx.LinxInterpretedData;
import com.hartwig.hmftools.datamodel.cohort.PercentileType;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleInterpretedData;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
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
    public abstract LinxInterpretedData linx();

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
    public abstract List<PeachGenotype> peach();

    @Nullable
    public abstract List<SignatureAllocation> sigAllocations();

    @NotNull
    public abstract Map<PercentileType, Evaluation> cohortEvaluations();

    @NotNull
    public abstract OrangePlots plots();

}
