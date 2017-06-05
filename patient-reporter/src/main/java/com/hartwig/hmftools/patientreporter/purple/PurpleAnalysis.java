package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleAnalysis {

    @NotNull
    public abstract FittedPurity fittedPurity();

    @NotNull
    public abstract List<PurpleCopyNumber> copyNumbers();

    @NotNull
    public List<VariantReport> enrich(@NotNull List<VariantReport> variants) {
        return VariantCopyNumberZipper.zip(fittedPurity(), variants, copyNumbers());
    }
}
