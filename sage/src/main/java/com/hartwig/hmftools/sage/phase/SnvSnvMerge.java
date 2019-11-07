package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class SnvSnvMerge implements Consumer<SageVariant> {

    private static final int BUFFER = 2;

    private final MnvFactory factory;
    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> deque = Lists.newLinkedList();

    SnvSnvMerge(@NotNull final Consumer<SageVariant> consumer, MnvFactory factory) {
        this.consumer = consumer;
        this.factory = factory;
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        flush(newEntry);
        if (newEntry.isPassing() && newEntry.localPhaseSet() > 0 && !newEntry.isIndel()) {

            for (int i = 0; i < deque.size(); i++) {
                final SageVariant oldEntry = deque.get(i);
                if (oldEntry.isPassing() && oldEntry.localPhaseSet() == newEntry.localPhaseSet() && !newEntry.isIndel()) {
                    SageVariant mnv = factory.createMNV(oldEntry, newEntry);
                    deque.add(i, mnv);
                    if (mnv.isPassing()) {
                        oldEntry.filters().add("merge");
                        newEntry.filters().add("merge");
                    }
                    i++;
                }
            }
        }

        deque.add(newEntry);
    }

    private void flush(GenomePosition position) {
        Iterator<SageVariant> iterator = deque.iterator();
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
        deque.forEach(consumer);
        deque.clear();
    }

}
