package com.hartwig.hmftools.ckb.datamodel;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntry {

    public abstract int profileId();

    @NotNull
    public abstract LocalDate createDate();

    @NotNull
    public abstract LocalDate updateDate();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<Variant> variants();

    @NotNull
    public abstract List<Evidence> evidences();

    @NotNull
    public abstract List<ClinicalTrial> clinicalTrials();
}