package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.Nullable;

public interface GermlineVariant extends VariantProxy
{
    @Nullable
    String pathogenicity();

    boolean pathogenic();
}
