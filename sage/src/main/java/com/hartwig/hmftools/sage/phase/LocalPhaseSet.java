package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class LocalPhaseSet extends BufferedPostProcessor {

    private int phase;

    LocalPhaseSet(int flankSize, @NotNull final Consumer<SageVariant> consumer) {
        super(flankSize + 25, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newEntry, @NotNull final Collection<SageVariant> buffer) {
        final ReadContext newReadContext = newEntry.readContext();
        for (final SageVariant oldEntry : buffer) {
            final ReadContext oldReadContext = oldEntry.readContext();
            int offset = offset(oldEntry.variant(), newEntry.variant());

            if (oldReadContext.phased(offset, newReadContext)) {
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
    }

    static int positionOffset(@NotNull final VariantHotspot left, @NotNull final VariantHotspot right) {
        long positionOffset = left.position() - right.position();
        return (int) (positionOffset);
    }

    static int offset(@NotNull final VariantHotspot left, @NotNull final VariantHotspot right) {

        long positionOffset = positionOffset(left, right);
        if (positionOffset == 0) {
            return 0;
        }

        return (int) (positionOffset + Math.max(0, left.ref().length() - left.alt().length()) - Math.max(0,
                left.alt().length() - left.ref().length()));
    }

}
