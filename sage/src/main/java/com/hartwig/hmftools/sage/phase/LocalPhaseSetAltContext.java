package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class LocalPhaseSetAltContext implements Consumer<AltContext> {

    private int phase;
    private final Consumer<AltContext> consumer;
    private final ArrayDeque<AltContext> deque = new ArrayDeque<>();

    public LocalPhaseSetAltContext(@NotNull final Consumer<AltContext> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final AltContext newEntry) {
        final ReadContext newReadContext = newEntry.primaryReadContext().readContext();;

        Iterator<AltContext> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final AltContext oldEntry = iterator.next();

            final ReadContext oldReadContext = oldEntry.primaryReadContext().readContext();
            long distance = newEntry.position() - oldEntry.position();
            int offset = offset(oldEntry, newEntry);

            if (distance > 30 || distance < 0) {
                iterator.remove();
                consumer.accept(oldEntry);
            } else if (oldReadContext.phased(offset, newReadContext)) {
                if (oldEntry.localPhaseSet() != 0) {
                    newEntry.localPhaseSet(oldEntry.localPhaseSet());
                } else if (newEntry.localPhaseSet() != 0) {
                    oldEntry.localPhaseSet(newEntry.localPhaseSet());
                } else {
                    phase++;
                    oldEntry.localPhaseSet(phase);
                    newEntry.localPhaseSet(phase);
                }
            }
        }

        deque.add(newEntry);
    }


    public void flush() {
        deque.forEach(consumer);
        deque.clear();
    }

    private static int offset(@NotNull final VariantHotspot left, @NotNull final VariantHotspot right) {

        long positionOffset = left.position() - right.position();
        if (positionOffset == 0) {
            return 0;
        }

        return (int) (positionOffset + Math.max(0, left.ref().length() - left.alt().length()) - Math.max(0,
                left.alt().length() - left.ref().length()));
    }

}
