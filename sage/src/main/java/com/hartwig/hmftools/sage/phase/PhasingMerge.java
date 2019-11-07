package com.hartwig.hmftools.sage.phase;

import java.util.function.Consumer;

import com.hartwig.hmftools.sage.variant.SageVariant;

public class PhasingMerge implements Consumer<SageVariant> {

    private final Consumer<SageVariant> consumer;

    public PhasingMerge(final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(final SageVariant entry) {

    }
}
