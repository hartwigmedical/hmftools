package com.hartwig.hmftools.ckb.datamodel.common.therapy;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.drug.Drug;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Therapy {

    public abstract int id();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract List<Drug> drugs();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<TherapyDescription> descriptions();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses();
}