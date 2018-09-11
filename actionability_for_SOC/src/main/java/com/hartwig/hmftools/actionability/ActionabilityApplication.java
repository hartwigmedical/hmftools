package com.hartwig.hmftools.actionability;

import org.apache.logging.log4j.LogManager;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityApplication {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityApplication.class);

    public static void main(final String... args)  {
        LOGGER.info("hello");

    }




}
