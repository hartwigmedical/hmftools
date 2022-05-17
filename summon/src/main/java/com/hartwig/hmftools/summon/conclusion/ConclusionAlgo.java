package com.hartwig.hmftools.summon.conclusion;

import java.io.IOException;

import com.hartwig.hmftools.summon.SummonData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionAlgo {

    @NotNull
    public static ActionabilityConclusion generateConclusion(@NotNull SummonData summonData) throws IOException {
        return ImmutableActionabilityConclusion.builder().conclusion(Strings.EMPTY).build();
    }
}
