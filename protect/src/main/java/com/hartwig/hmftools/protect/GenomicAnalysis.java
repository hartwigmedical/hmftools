package com.hartwig.hmftools.protect;

import java.util.List;

import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertion;

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
    public abstract List<EvidenceItem> tumorSpecificEvidence();

    @NotNull
    public abstract List<ClinicalTrial> clinicalTrials();

    @NotNull
    public abstract List<EvidenceItem> offLabelEvidence();

    @NotNull
    public abstract List<ReportableVariant> reportableVariants();

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
    public abstract List<LinxFusion> geneFusions();

    @NotNull
    public abstract List<ReportableGeneDisruption> geneDisruptions();

    @NotNull
    public abstract List<ReportableHomozygousDisruption> homozygousDisruptions();

    @Nullable
    public abstract List<ViralInsertion> viralInsertions();
}
