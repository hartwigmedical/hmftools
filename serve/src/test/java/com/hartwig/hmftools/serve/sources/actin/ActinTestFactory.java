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
        return ImmutableActinEntry.builder().trial(Strings.EMPTY).rule(ActinRule.ACTIVATION_OF_GENE_X).gene(Strings.EMPTY).build();
    }
}
