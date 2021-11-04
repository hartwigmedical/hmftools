package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GenomicAnalysis {

    public abstract double impliedPurity();

    public abstract boolean hasReliablePurity();

    public abstract boolean hasReliableQuality();

    public abstract double averageTumorPloidy();

    @NotNull
    public abstract List<ProtectEvidence> tumorSpecificEvidence();

    @NotNull
    public abstract List<ProtectEvidence> clinicalTrials();

    @NotNull
    public abstract List<ProtectEvidence> offLabelEvidence();

    @NotNull
    public abstract List<ReportableVariant> reportableVariants();

    @NotNull
    public abstract Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant();

    public abstract double microsatelliteIndelsPerMb();

    @NotNull
    public abstract MicrosatelliteStatus microsatelliteStatus();

    public abstract int tumorMutationalLoad();

    @NotNull
    public abstract TumorMutationalStatus tumorMutationalLoadStatus();

    public abstract double tumorMutationalBurden();

    public abstract double chordHrdValue();

    @NotNull
    public abstract ChordStatus chordHrdStatus();

    @NotNull
    public abstract List<ReportableGainLoss> gainsAndLosses();

    @NotNull
    public abstract List<CnPerChromosomeArmData> cnPerChromosome();

    @NotNull
    public abstract List<LinxFusion> geneFusions();

    @NotNull
    public abstract List<ReportableGeneDisruption> geneDisruptions();

    @NotNull
    public abstract List<ReportableHomozygousDisruption> homozygousDisruptions();

    @NotNull
    public abstract List<AnnotatedVirus> reportableViruses();
}
