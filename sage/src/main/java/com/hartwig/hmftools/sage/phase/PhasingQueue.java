package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class PhasingQueue implements Consumer<SageVariant> {

    private int phase;
    private final Consumer<SageVariant> consumer;
    private final ArrayDeque<SageVariant> deque = new ArrayDeque<>();

    public PhasingQueue(@NotNull final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final SageVariant entry) {
        final AltContext newAltContext = entry.primaryTumor();
        final ReadContext newReadContext = newAltContext.primaryReadContext().readContext();

        Iterator<SageVariant> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final SageVariant oldEntry = iterator.next();

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
