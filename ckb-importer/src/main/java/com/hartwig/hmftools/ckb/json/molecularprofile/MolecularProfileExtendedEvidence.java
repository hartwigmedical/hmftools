package com.hartwig.hmftools.ckb.json.molecularprofile;

import java.util.List;

import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileExtendedEvidence {

    public abstract int totalCount();

    @NotNull
    public abstract List<EvidenceInfo> evidence();
}
