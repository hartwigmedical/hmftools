package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.count.EvictingArray;

import org.jetbrains.annotations.NotNull;

public class TumorRefContextCandidates implements RefContextCandidates {

    private final String sample;
    private final EvictingArray<RefContext> rollingCandidates;
    private final List<RefContext> savedCandidates = Lists.newArrayList();

    public TumorRefContextCandidates(@NotNull final String sample) {
        this.sample = sample;
        final Consumer<RefContext> evictionHandler = (refContext) -> {
            if (!refContext.isAltsEmpty()) {
                savedCandidates.add(refContext);
            }
        };

        this.rollingCandidates = new EvictingArray<>(256, evictionHandler);
    }

    @NotNull
    public RefContext refContext(@NotNull final String chromosome, final long position) {
        return rollingCandidates.computeIfAbsent(position, aLong -> new RefContext(sample, chromosome, position));
    }

    @NotNull
    public List<RefContext> refContexts() {
        rollingCandidates.evictAll();
        Collections.sort(savedCandidates);
        return savedCandidates;
    }

}
