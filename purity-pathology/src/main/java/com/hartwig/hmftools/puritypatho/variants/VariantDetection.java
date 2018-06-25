package com.hartwig.hmftools.puritypatho.variants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public class VariantDetection {
    private static final Logger LOGGER = LogManager.getLogger(VariantDetection.class);

}
