package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberAnalysis {

    public abstract double purity();

    public abstract boolean hasReliablePurity();

    public abstract boolean hasReliableQuality();

    public abstract double ploidy();

    @NotNull
    public abstract List<GeneCopyNumber> exomeGeneCopyNumbers();

    @NotNull
    public abstract List<ReportableGainLoss> reportableGainsAndLosses();

    @NotNull
    public abstract List<EvidenceItem> evidenceItems();

    @NotNull
    public abstract PurpleSignatures purpleSignatures();

}
