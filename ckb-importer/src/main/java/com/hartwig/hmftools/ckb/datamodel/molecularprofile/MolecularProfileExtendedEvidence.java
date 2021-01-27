package com.hartwig.hmftools.ckb.datamodel.molecularprofile;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileExtendedEvidence {

    @NotNull
    public abstract String totalCount();

    @NotNull
    public abstract List<EvidenceInfo> evidence();
}
