package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class NormalRefContextCandidates implements RefContextCandidates {

    private final String sample;
    private final List<RefContext> refContexts;
    private final GenomePositionSelector<RefContext> refContextSelector;

    public NormalRefContextCandidates(@NotNull final String sample) {
        this.sample = sample;
        this.refContexts = Lists.newArrayList();
        this.refContextSelector = GenomePositionSelectorFactory.create(refContexts);
    }

    @NotNull
    public RefContext add(@NotNull final String chromosome, final long position) {
        if (!refContexts.isEmpty() && refContexts.get(refContexts.size() - 1).position() > position) {
            throw new IllegalArgumentException("Can only add sorted ref contexts");
        }

        Optional<RefContext> optionalRefContext = refContextSelector.select(GenomePositions.create(chromosome, position));
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
        return refContextSelector.select(GenomePositions.create(chromosome, position)).orElse(null);
    }

    @NotNull
    @Override
    public List<RefContext> refContexts() {
        return refContexts;
    }
}
