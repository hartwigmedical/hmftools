package com.hartwig.hmftools.serve.refgenome.liftover;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface LiftOverAlgo {

    @Nullable
    LiftOverResult liftOver(@NotNull String chromosome, int position);
}
