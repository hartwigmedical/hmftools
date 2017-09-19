package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantAnalysis {

    @NotNull
    public abstract List<SomaticVariant> allVariants();

    @NotNull
    public abstract List<SomaticVariant> passedVariants();

    @NotNull
    public abstract List<SomaticVariant> missenseVariants();

    @NotNull
    public abstract List<SomaticVariant> consequentialVariants();

    @NotNull
    public abstract List<SomaticVariant> potentialConsequentialMNVs();

    @NotNull
    public abstract List<VariantReport> findings();

    public int mutationalLoad() {
        return missenseVariants().size();
    }
}
