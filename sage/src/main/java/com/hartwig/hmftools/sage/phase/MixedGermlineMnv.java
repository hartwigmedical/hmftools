package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class MixedGermlineMnv implements Consumer<SageVariant> {

    private static final int BUFFER = 10;

    private final ArrayDeque<SageVariant> buffer = new ArrayDeque<>();
    private final Consumer<SageVariant> consumer;

    public MixedGermlineMnv(final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(final SageVariant newVariant) {
        flush(newVariant);
        buffer.add(newVariant);

        if (!newVariant.isIndel()) {
            boolean newIsPassingMnv = isPassingMnv(newVariant);
            boolean newIsGermlineSnv = isGermlineFilteredSnv(newVariant);
            if (newIsPassingMnv || newIsGermlineSnv) {
                for (final SageVariant other : buffer) {
                    if (newIsPassingMnv && isGermlineFilteredSnv(other)) {
                        if (DedupMnv.longerContainsShorter(other, newVariant)) {
                            newVariant.filters().add(SoftFilter.MIXED_GERMLINE_SOMATIC_MNV.toString());
                        }
                    } else if (newIsGermlineSnv && isPassingMnv(other)) {
                        if (DedupMnv.longerContainsShorter(newVariant, other)) {
                            other.filters().add(SoftFilter.MIXED_GERMLINE_SOMATIC_MNV.toString());
                        }
                    }
                }
            }
        }
    }

    private static boolean isGermlineFilteredSnv(final SageVariant newVariant) {
        return !newVariant.isPassing() && SoftFilter.isGermlineAndNotTumorFiltered(newVariant.filters()) && newVariant.ref().length() == 1;
    }

    private static boolean isPassingMnv(final SageVariant newVariant) {
        return newVariant.isPassing() && newVariant.ref().length() > 1;
    }

    public void flush() {
        buffer.forEach(consumer);
        buffer.clear();
    }

    private void flush(@NotNull final GenomePosition position) {
        final Iterator<SageVariant> iterator = buffer.iterator();
        while (iterator.hasNext()) {
            final SageVariant entry = iterator.next();
            long entryEnd = entry.position() + entry.ref().length() - 1;
            if (!entry.chromosome().equals(position.chromosome()) || entryEnd < position.position() - BUFFER) {
                iterator.remove();
                consumer.accept(entry);
            } else {
                return;
            }
        }
    }

}
