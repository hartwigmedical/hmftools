package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberAnalysis {

    public abstract boolean hasReliablePurityFit();

    public abstract double purity();

    public abstract double ploidy();

    @NotNull
    public abstract Gender gender();

    @NotNull
    public abstract List<GeneCopyNumber> exomeGeneCopyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> reportableGeneCopyNumbers();

    @NotNull
    public abstract List<EvidenceItem> evidenceItems();

}
