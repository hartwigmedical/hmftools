package com.hartwig.hmftools.sage.count;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.ReadContext;

import org.jetbrains.annotations.NotNull;

public class PhasingQueue implements Consumer<AltContext> {

    private int phase;
    private final Consumer<AltContext> consumer;
    private final ArrayDeque<AltContext> deque = new ArrayDeque<>();

    public PhasingQueue(@NotNull final Consumer<AltContext> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(final AltContext newAltContext) {
        final ReadContext newReadContext = newAltContext.primaryReadContext().readContext();

        Iterator<AltContext> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final AltContext oldAltContext = iterator.next();
            final ReadContext oldReadContext = oldAltContext.primaryReadContext().readContext();
            long distance = newAltContext.position() - oldAltContext.position();

            if (distance > 30 || distance < 0) {
                iterator.remove();
                consumer.accept(oldAltContext);
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

        deque.add(newAltContext);
    }

    public void flush() {
        deque.forEach(consumer);
        deque.clear();
    }

}
