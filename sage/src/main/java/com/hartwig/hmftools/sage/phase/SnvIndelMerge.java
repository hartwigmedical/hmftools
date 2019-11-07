package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class SnvIndelMerge implements Consumer<SageVariant> {

    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> list = Lists.newLinkedList();

    SnvIndelMerge(@NotNull final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        flush(newEntry);
        if (newEntry.isPassing() && newEntry.localPhaseSet() > 0) {

            for (int i = 0; i < list.size(); i++) {
                final SageVariant oldEntry = list.get(i);
                if (oldEntry.isPassing() && oldEntry.localPhaseSet() == newEntry.localPhaseSet()
                        && oldEntry.position() == newEntry.position()) {

                    if (newEntry.isIndel() && !oldEntry.isIndel()) {
                        removeSnvInIndel(oldEntry, newEntry);
                    } else if (!newEntry.isIndel() && oldEntry.isIndel()) {
                        removeSnvInIndel(newEntry, oldEntry);
                    }

                    i++;
                }
            }
        }

        list.add(newEntry);
    }

    private void removeSnvInIndel(@NotNull final SageVariant snv, @NotNull final SageVariant indel) {
        if (!snv.isSynthetic() && snv.normal().alt().equals(indel.normal().alt().substring(0, 1))) {
            snv.filters().add("merge");
        }
    }

    private void flush(GenomePosition position) {
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

}
