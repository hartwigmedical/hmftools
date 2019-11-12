package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.SageVCF;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class DedupIndel implements Consumer<SageVariant> {

    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> list = Lists.newArrayList();

    public DedupIndel(final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(final SageVariant variant) {
        flush(variant);
        if (isPassingPhasedIndel(variant)) {
            for (final SageVariant other : list) {
                if (isPassingPhasedIndel(other) && variant.localPhaseSet() == other.localPhaseSet()) {

                    if (variant.isDelete() && other.isDelete()) {
                        processDel(variant, other);
                    }

                    if (variant.isInsert() && other.isInsert()) {
                        processIns(variant, other);
                    }
                }
            }
        }

        list.add(variant);

    }

    private void processDel(@NotNull final SageVariant left, @NotNull final SageVariant right) {
        if (!left.normal().alt().equals(right.normal().alt())) {
            return;
        }

        final SageVariant shorter;
        final SageVariant longer;
        if (left.normal().ref().length() < right.normal().ref().length()) {
            shorter = left;
            longer = right;
        } else {
            shorter = right;
            longer = left;
        }

        if (longer.normal().ref().substring(0, shorter.normal().ref().length()).equals(shorter.normal().ref())) {
            longer.filters().add(SageVCF.DEDUP_FILTER);
        }
    }

    private void processIns(@NotNull final SageVariant left, @NotNull final SageVariant right) {
        if (!left.normal().ref().equals(right.normal().ref())) {
            return;
        }

        final SageVariant shorter;
        final SageVariant longer;
        if (left.normal().alt().length() < right.normal().alt().length()) {
            shorter = left;
            longer = right;
        } else {
            shorter = right;
            longer = left;
        }

        if (longer.normal().alt().substring(0, shorter.normal().alt().length()).equals(shorter.normal().alt())) {
            longer.filters().add(SageVCF.DEDUP_FILTER);
        }

    }

    private void flush(@NotNull final GenomePosition position) {
        final Iterator<SageVariant> iterator = list.iterator();
        while (iterator.hasNext()) {
            final SageVariant entry = iterator.next();
            if (!entry.chromosome().equals(position.chromosome()) || entry.position() < position.position()) {
                iterator.remove();
                consumer.accept(entry);
            } else {
                return;
            }
        }
    }

    public void flush() {
        list.forEach(consumer);
        list.clear();
    }

    private boolean isPassingPhasedIndel(@NotNull final SageVariant newEntry) {
        return newEntry.isPassing() && newEntry.localPhaseSet() > 0 && newEntry.isIndel();
    }

}
