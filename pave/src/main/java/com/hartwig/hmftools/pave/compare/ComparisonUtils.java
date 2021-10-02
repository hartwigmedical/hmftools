package com.hartwig.hmftools.pave.compare;

import java.util.List;

import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.pave.VariantTransImpact;

public final class ComparisonUtils
{
    public static SnpEffAnnotation findMatchingAnnotation(final VariantTransImpact transImpact, final List<SnpEffAnnotation> annotations)
    {
        return annotations.stream()
                .filter(x -> x.featureID() != null && x.featureID().equals(transImpact.TransData.TransName))
                .findFirst().orElse(null);
    }

}
