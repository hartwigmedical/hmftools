package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.SageVCF;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class SnvSnvMerge implements Consumer<SageVariant> {

    private static final int BUFFER = 2;

    private final MnvFactory factory;
    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> list = Lists.newLinkedList();

    SnvSnvMerge(@NotNull final Consumer<SageVariant> consumer, MnvFactory factory) {
        this.consumer = consumer;
        this.factory = factory;
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        flush(newEntry);
        if (isEligibleWithAny(newEntry)) {

            for (int i = 0; i < list.size(); i++) {
                final SageVariant oldEntry = list.get(i);
                if (isEligibleWithLps(newEntry.localPhaseSet(), oldEntry)) {
                    SageVariant mnv = factory.createMNV(oldEntry, newEntry);
                    list.add(i, mnv);
                    if (mnv.isPassing()) {
                        oldEntry.filters().add(SageVCF.MERGE_FILTER);
                        if (oldEntry.isSynthetic()) {
                            list.remove(i + 1);
                        }
                        i--;
                        newEntry.filters().add(SageVCF.MERGE_FILTER);
                    }
                    i++;
                }
            }
        }

        list.add(newEntry);
    }

    private void flush(@NotNull final GenomePosition position) {
        final Iterator<SageVariant> iterator = list.iterator();
        while (iterator.hasNext()) {
            final SageVariant entry = iterator.next();
            long entryEnd = entry.position() + entry.normal().ref().length() - 1;
            if (!entry.chromosome().equals(position.chromosome()) || entryEnd < position.position() - BUFFER) {
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

    private boolean isEligibleWithAny(@NotNull final SageVariant entry) {
        return entry.isPassing() && entry.localPhaseSet() > 0 && !entry.isIndel();
    }

    private boolean isEligibleWithLps(int lps, @NotNull final SageVariant entry) {
        return isEligibleWithAny(entry) && entry.localPhaseSet() == lps;
    }

}
