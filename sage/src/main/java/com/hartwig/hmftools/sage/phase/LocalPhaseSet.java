package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class LocalPhaseSet implements Consumer<SageVariant> {

    private int phase;
    private final boolean germline;
    private final Consumer<SageVariant> consumer;
    private final ArrayDeque<SageVariant> deque = new ArrayDeque<>();

    LocalPhaseSet(boolean germline, @NotNull final Consumer<SageVariant> consumer) {
        this.germline = germline;
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        final AltContext newAltContext = altContext(newEntry);
        final ReadContext newReadContext = newAltContext.primaryReadContext().readContext();;

        Iterator<SageVariant> iterator = deque.iterator();
        while (iterator.hasNext()) {
            final SageVariant oldEntry = iterator.next();

            final AltContext oldAltContext = altContext(oldEntry);
            final ReadContext oldReadContext = oldAltContext.primaryReadContext().readContext();
            long distance = newAltContext.position() - oldAltContext.position();
            int offset = offset(oldEntry.normal(), newEntry.normal());

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

    private AltContext altContext(@NotNull final SageVariant newEntry) {
        return germline ? newEntry.normal() : newEntry.primaryTumor();
    }

    public void flush() {
        deque.forEach(consumer);
        deque.clear();
    }

    static int offset(@NotNull final VariantHotspot left, @NotNull final VariantHotspot right) {

        long positionOffset = left.position() - right.position();
        if (positionOffset == 0) {
            return 0;
        }

        return (int) (positionOffset + Math.max(0, left.ref().length() - left.alt().length()) - Math.max(0,
                left.alt().length() - left.ref().length()));
    }

}
