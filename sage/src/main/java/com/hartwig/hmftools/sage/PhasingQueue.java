package com.hartwig.hmftools.sage;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class PhasingQueue implements Consumer<SageEntry> {

    private int phase;
    private final Consumer<SageEntry> consumer;
    private final ArrayDeque<SageEntry> deque = new ArrayDeque<>();

    PhasingQueue(@NotNull final Consumer<SageEntry> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final SageEntry entry) {
        final AltContext newAltContext = entry.primaryTumor();
        final ReadContext newReadContext = newAltContext.primaryReadContext().readContext();

        Iterator<SageEntry> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final SageEntry oldEntry = iterator.next();

            final AltContext oldAltContext = oldEntry.primaryTumor();
            final ReadContext oldReadContext = oldAltContext.primaryReadContext().readContext();
            long distance = newAltContext.position() - oldAltContext.position();

            if (distance > 30 || distance < 0) {
                iterator.remove();
                consumer.accept(oldEntry);
            } else if (oldReadContext.phased(newReadContext)) {
                if (oldAltContext.phase() != 0) {
                    newAltContext.phase(oldAltContext.phase());
                } else if (newAltContext.phase() != 0) {
                    oldAltContext.phase(newAltContext.phase());
                } else {
                    phase++;
                    oldAltContext.phase(phase);
                    newAltContext.phase(phase);
                }
            }
        }

        deque.add(entry);
    }

    public void flush() {
        deque.forEach(consumer);
        deque.clear();
    }

}
