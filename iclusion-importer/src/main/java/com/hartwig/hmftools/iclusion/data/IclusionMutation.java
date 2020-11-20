package com.hartwig.hmftools.iclusion.data;

import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.iclusion.classification.MutationTypeClassifier;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionMutation {

    @NotNull
    @Value.Derived
    public MutationType type() {
        return MutationTypeClassifier.classify(this);
    }

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String name();

    public abstract boolean negation();
}
