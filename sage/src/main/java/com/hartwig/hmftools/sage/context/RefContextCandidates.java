package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.count.EvictingLinkedMap;

import org.jetbrains.annotations.NotNull;

public class RefContextCandidates {

    public static RefContextCandidates tumorCandidates() {
        return new RefContextCandidates();
    }

    public static RefContextCandidates normalCandidates(@NotNull final Set<Long> hotspots) {
        return new RefContextCandidates(hotspots);
    }

    private final EvictingLinkedMap<Long, RefContext> rollingCandidates;
    private final List<RefContext> savedCandidates = Lists.newArrayList();

    private RefContextCandidates(@NotNull final Set<Long> hotspots) {

        final BiConsumer<Long, RefContext> evictionHandler = (position, refContext) -> {
            if (hotspots.contains(position)) {
                savedCandidates.add(refContext);
            }
        };

        this.rollingCandidates = new EvictingLinkedMap<>(evictionHandler);
    }

    private RefContextCandidates() {
        final BiConsumer<Long, RefContext> evictionHandler = (position, refContext) -> {
            if (!refContext.isAltsEmpty()) {
                savedCandidates.add(refContext);
            }
        };

        this.rollingCandidates = new EvictingLinkedMap<>(evictionHandler);
    }

    @NotNull
    public RefContext refContext(@NotNull final String sample, @NotNull final String chromosome, final long position,
            @NotNull final String ref) {
        return rollingCandidates.computeIfAbsent(position, aLong -> new RefContext(sample, chromosome, position, ref));
    }

    @NotNull
    public List<RefContext> refContexts() {
        rollingCandidates.evictAll();
        Collections.sort(savedCandidates);
        return savedCandidates;
    }

    @NotNull
    public List<AltContext> altContexts() {
        return refContexts().stream().flatMap(x -> x.alts().stream()).collect(Collectors.toList());
    }

}
