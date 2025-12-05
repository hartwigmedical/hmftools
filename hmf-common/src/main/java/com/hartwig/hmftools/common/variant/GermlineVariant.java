package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.Nullable;

public interface GermlineVariant extends VariantDelegate
{
    @Nullable
    String pathogenicity();

    boolean pathogenic();
}
