package com.hartwig.hmftools.ckb.molecularprofile;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileExtendedEvidence {

    @NotNull
    public abstract String totalCount();

    @NotNull
    public abstract List<MolecularProfileExtendedEvidenceList> extendedEvidenceList();
}
