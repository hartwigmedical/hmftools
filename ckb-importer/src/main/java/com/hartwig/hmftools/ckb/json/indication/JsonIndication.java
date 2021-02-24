package com.hartwig.hmftools.ckb.json.indication;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonIndication implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String source();

    @Nullable
    public abstract String definition();

    @Nullable
    public abstract String currentPreferredTerm();

    @Nullable
    public abstract LocalDate lastUpdateDateFromDO();

    @NotNull
    public abstract List<String> altIds();

    @NotNull
    public abstract String termId();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrials();
}
