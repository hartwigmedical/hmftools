package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActinTestFactory {

    private ActinTestFactory() {
    }

    @NotNull
    public static ActinEntry createTestEntry() {
        return ImmutableActinEntry.builder()
                .trial(Strings.EMPTY)
                .rule(ActinRule.ACTIVATION_OR_AMPLIFICATION_OF_GENE_X)
                .gene(Strings.EMPTY)
                .mutation(Strings.EMPTY)
                .build();
    }

    @NotNull
    public static ActinEntry createTestEntryWithData(@NotNull ActinRule actinRule, @NotNull String gene, @NotNull String mutation) {
        return ImmutableActinEntry.builder()
                .trial("A")
                .rule(actinRule)
                .gene(gene)
                .mutation(mutation)
                .build();
    }
}