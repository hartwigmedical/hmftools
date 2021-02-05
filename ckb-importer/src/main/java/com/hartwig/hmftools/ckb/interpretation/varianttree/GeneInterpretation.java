package com.hartwig.hmftools.ckb.interpretation.varianttree;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneInterpretation {

    @NotNull
    public abstract Gene gene();

    @NotNull
    public abstract List<Reference> references();
}
