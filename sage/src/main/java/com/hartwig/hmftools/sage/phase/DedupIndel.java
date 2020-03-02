package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;

public class DedupIndel implements Consumer<SageVariant> {

    private final Consumer<SageVariant> consumer;
    private final ArrayDeque<SageVariant> list =  new ArrayDeque<>();

    DedupIndel(final Consumer<SageVariant> consumer) {
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
        if (!left.alt().equals(right.alt())) {
            return;
        }

        final SageVariant shorter;
        final SageVariant longer;
        if (left.ref().length() < right.ref().length()) {
            shorter = left;
            longer = right;
        } else {
            shorter = right;
            longer = left;
        }

        if (longer.ref().substring(0, shorter.ref().length()).equals(shorter.ref())) {
            longer.filters().add(SageVCF.DEDUP_FILTER);
        }
    }

    private void processIns(@NotNull final SageVariant left, @NotNull final SageVariant right) {
        if (!left.ref().equals(right.ref())) {
            return;
        }

        final SageVariant shorter;
        final SageVariant longer;
        if (left.alt().length() < right.alt().length()) {
            shorter = left;
            longer = right;
        } else {
            shorter = right;
            longer = left;
        }

        if (longer.alt().substring(0, shorter.alt().length()).equals(shorter.alt())) {
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
