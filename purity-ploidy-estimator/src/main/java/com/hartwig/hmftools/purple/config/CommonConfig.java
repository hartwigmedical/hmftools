package com.hartwig.hmftools.purple.config;

import java.util.Map;

import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.refgenome.RefGenome;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CommonConfig {

    int WINDOW_SIZE = 1000;

    @NotNull
    String refSample();

    @NotNull
    String tumorSample();

    @NotNull
    String outputDirectory();

    @NotNull
    String cobaltDirectory();

    @NotNull
    String amberDirectory();

    @NotNull
    String gcProfile();

    default int windowSize() {
        return WINDOW_SIZE;
    }
}
