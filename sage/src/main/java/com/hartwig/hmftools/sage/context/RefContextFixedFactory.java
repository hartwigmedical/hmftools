package com.hartwig.hmftools.sage.context;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RefContextFixedFactory implements RefContextFactory {

    private final String sample;
    private final List<RefContextFixed> refContexts;
    private final PositionSelector<RefContextFixed> refPositionSelector;

    RefContextFixedFactory(@NotNull final String sample) {
        this.sample = sample;
        this.refContexts = Lists.newArrayList();
        this.refPositionSelector = new PositionSelector<>(refContexts);
    }

    public void create(@NotNull final String chromosome, final long position, final String ref, final String alt,
            final ReadContext readContext) {
        final RefContextFixed refContext = add(chromosome, position);
        refContext.altContext(ref, alt, readContext);
    }

    @NotNull
    private RefContextFixed add(@NotNull final String chromosome, final long position) {
        if (!refContexts.isEmpty() && refContexts.get(refContexts.size() - 1).position() > position) {
            throw new IllegalArgumentException("Can only add sorted ref contexts");
        }

        Optional<RefContextFixed> optionalRefContext = refPositionSelector.select(position);
        if (optionalRefContext.isPresent()) {
            return optionalRefContext.get();
        }

        RefContextFixed refContext = new RefContextFixed(sample, chromosome, position);
        refContexts.add(refContext);

        return refContext;
    }

    @Nullable
    @Override
    public RefContext refContext(@NotNull final String chromosome, final long position) {
        return refPositionSelector.select(position).orElse(null);
    }

    @NotNull
    public List<RefContext> refContexts() {
        return new ArrayList<>(refContexts);
    }
}
