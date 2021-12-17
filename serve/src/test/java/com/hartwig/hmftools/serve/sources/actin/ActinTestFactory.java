package com.hartwig.hmftools.serve.sources.actin;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.jetbrains.annotations.NotNull;

public class ActinTestFactory {

    private ActinTestFactory() {
    }

    @NotNull
    public static ActinEntry createEntry() {
        return ImmutableActinEntry.builder()
                .trial("trial 1")
                .rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y)
                .parameters(Lists.newArrayList("B", "remove"))
                .build();
    }
}
