package com.hartwig.hmftools.purple.config;

import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CircosConfig {
    public abstract Optional<String> circosBinary();

    public abstract String plotDirectory();

    public abstract String circosDirectory();
}
