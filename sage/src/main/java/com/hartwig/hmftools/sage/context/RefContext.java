package com.hartwig.hmftools.sage.context;

import java.util.Collection;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public interface RefContext extends GenomePosition {

    @NotNull
    String sample();

    int rawDepth();

    void refRead();

    void altReadFixed(@NotNull final String ref, @NotNull final String alt, int baseQuality);

    void altReadCandidate(@NotNull final String ref, @NotNull final String alt, int baseQuality, @NotNull final ReadContext interimReadContext);

    @NotNull
    Collection<AltContext> alts();

}
