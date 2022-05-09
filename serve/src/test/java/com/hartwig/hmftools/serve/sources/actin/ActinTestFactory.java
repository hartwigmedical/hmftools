package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActinTestFactory {

    private ActinTestFactory() {
    }

    @NotNull
    public static ImmutableActinEntry.Builder builder() {
        return ImmutableActinEntry.builder()
                .trial(Strings.EMPTY)
                .rule(ActinRule.ACTIVATION_OR_AMPLIFICATION_OF_GENE_X)
                .gene(Strings.EMPTY)
                .mutation(Strings.EMPTY)
                .isUsedAsInclusion(false);
    }

    @NotNull
    public static ActinEntry create(@NotNull ActinRule actinRule, @Nullable String gene, @Nullable String mutation) {
        return builder().trial("A").rule(actinRule).gene(gene).mutation(mutation).isUsedAsInclusion(true).build();
    }
}