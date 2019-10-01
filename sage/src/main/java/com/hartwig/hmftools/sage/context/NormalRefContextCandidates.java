package com.hartwig.hmftools.sage.context;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class NormalRefContextCandidates implements RefContextCandidates {

    private final List<RefContext> refContexts;
    private final GenomePositionSelector<RefContext> refContextSelector;

    public NormalRefContextCandidates() {
        this.refContexts = Lists.newArrayList();
        this.refContextSelector = GenomePositionSelectorFactory.create(refContexts);
    }

    public void addRefContext(@NotNull final RefContext refContext) {
        refContexts.add(refContext);
    }

    @Nullable
    @Override
    public RefContext refContext(@NotNull final String chromosome, final long position, @NotNull final String ref) {
        return refContextSelector.select(GenomePositions.create(chromosome, position)).orElse(null);
    }

    @NotNull
    @Override
    public List<RefContext> refContexts() {
        return refContexts;
    }
}
