package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class MnvMerge implements Consumer<SageVariant> {

    private static final int BUFFER = 2;

    private int phase;
    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> deque = Lists.newLinkedList();

    MnvMerge(@NotNull final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        flush(newEntry.position() - BUFFER);
        if (newEntry.isPassing() && newEntry.localPhaseSet() > 0 && !newEntry.isIndel()) {

            for (int i = 0; i < deque.size(); i++) {
                final SageVariant oldEntry = deque.get(i);
                if (oldEntry.isPassing() && oldEntry.localPhaseSet() == newEntry.localPhaseSet() && !newEntry.isIndel()) {
//                    System.out.println("GHERE");

                }
            }
        }

        deque.add(newEntry);
    }

    private void flush(long position) {
        Iterator<SageVariant> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final SageVariant entry = iterator.next();
            if (entry.position() < position) {
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
