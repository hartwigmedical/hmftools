package com.hartwig.hmftools.ckb.datamodel.therapy;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Therapy {

    public abstract int id();

    @NotNull
    public abstract LocalDate createDate();

    @Nullable
    public abstract LocalDate updateDate();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract List<Drug> drugs();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<GlobalApprovalStatus> globalApprovalStatuses();

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<Reference> references();
}