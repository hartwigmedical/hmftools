package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.count.EvictingArray;

import org.jetbrains.annotations.NotNull;

public class RefContextFactory {

    private final String sample;
    private final EvictingArray<RefContext> rollingCandidates;
    private final List<AltContext> savedCandidates = Lists.newArrayList();

    public RefContextFactory(@NotNull final String sample) {
        this.sample = sample;
        final Consumer<RefContext> evictionHandler = (refContext) -> {
            refContext.alts()
                    .stream()
                    .filter(this::refPredicate)
                    .filter(x -> x.primaryReadContext().readContext().isComplete())
                    .forEach(savedCandidates::add);
        };

        this.rollingCandidates = new EvictingArray<>(256, evictionHandler);
    }

    @NotNull
    public RefContext refContext(@NotNull final String chromosome, final long position) {
        return rollingCandidates.computeIfAbsent(position, aLong -> new RefContext(sample, chromosome, position));
    }

    @NotNull
    public List<AltContext> altContexts() {
        rollingCandidates.evictAll();
        Collections.sort(savedCandidates);
        return savedCandidates;
    }

    private boolean refPredicate(@NotNull final AltContext altContext) {
        for (int i = 0; i < altContext.ref().length(); i++) {
            char base = altContext.ref().charAt(i);
            if (base != 'G' && base != 'A' && base != 'T' && base != 'C') {
                return false;
            }
        }

        return true;
    }

}
