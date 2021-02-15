package com.hartwig.hmftools.ckb.json.treatmentapproach;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonTreatmentApproach implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String profileName();

    @Nullable
    public abstract DrugClassInfo drugClass();

    @Nullable
    public abstract TherapyInfo therapy();

    @NotNull
    public abstract List<ReferenceInfo> references();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();
}
