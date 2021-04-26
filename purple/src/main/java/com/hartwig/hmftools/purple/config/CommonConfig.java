package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CommonConfig {

    @NotNull
    String refSample();

    @NotNull
    String tumorSample();

    @NotNull
    String outputDirectory();

    @NotNull
    String sampleDirectory();

    @NotNull
    String amberDirectory();

    @NotNull
    String cobaltDirectory();

    @NotNull
    String gcProfile();

    @NotNull
    String version();

    boolean tumorOnly();

    default int windowSize() {
        return WINDOW_SIZE;
    }
}
