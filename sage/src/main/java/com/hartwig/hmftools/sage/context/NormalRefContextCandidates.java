package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class NormalRefContextCandidates implements RefContextCandidates {

    private final String sample;
    private final List<RefContext> refContexts;
    private final PositionSelector<RefContext> refPositionSelector;

    public NormalRefContextCandidates(@NotNull final String sample) {
        this.sample = sample;
        this.refContexts = Lists.newArrayList();
        this.refPositionSelector = new PositionSelector<>(refContexts);
    }

    @NotNull
    public RefContext add(@NotNull final String chromosome, final long position) {
        if (!refContexts.isEmpty() && refContexts.get(refContexts.size() - 1).position() > position) {
            throw new IllegalArgumentException("Can only add sorted ref contexts");
        }

        Optional<RefContext> optionalRefContext = refPositionSelector.select(position);
        if (optionalRefContext.isPresent()) {
            return optionalRefContext.get();
        }

        RefContext refContext = new RefContext(sample, chromosome, position);
        refContexts.add(refContext);

        return refContext;
    }

    @Nullable
    @Override
    public RefContext refContext(@NotNull final String chromosome, final long position) {
        return refPositionSelector.select(position).orElse(null);
    }

    @NotNull
    @Override
    public List<RefContext> refContexts() {
        return refContexts;
    }
}
