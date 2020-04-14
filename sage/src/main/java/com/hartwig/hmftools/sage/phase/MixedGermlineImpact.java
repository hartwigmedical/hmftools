package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class MixedGermlineImpact implements Consumer<SageVariant> {

    private static final int BUFFER = 10;
    private final ArrayDeque<SageVariant> buffer = new ArrayDeque<>();
    private final Consumer<SageVariant> consumer;

    private int sharedImpact;

    public MixedGermlineImpact(final Consumer<SageVariant> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(final SageVariant newVariant) {
        flush(newVariant);

        int lps = newVariant.localPhaseSet();

        if (!newVariant.isIndel() && lps > 0) {

            boolean newVariantIsSnv = isPassingSnv(newVariant);
            boolean newVariantIsMixedGermlineMnv = isMixedGermlineMnv(newVariant);

            if (newVariantIsSnv || newVariantIsMixedGermlineMnv) {
                for (SageVariant oldVariant : buffer) {
                    if (oldVariant.localPhaseSet() == lps) {
                        if (newVariantIsSnv && isMixedGermlineMnv(oldVariant)) {
                            process(oldVariant, newVariant);
                        } else if (newVariantIsMixedGermlineMnv && isPassingSnv(oldVariant)) {
                            process(newVariant, oldVariant);
                        }
                    }
                }
            }

        }

        buffer.add(newVariant);
    }

    private void process(SageVariant mnv, SageVariant snv) {
        if (DedupMnv.longerContainsShorter(snv, mnv)) {

            if (mnv.mixedGermlineImpact() == 0) {
                mnv.mixedGermlineImpact(++sharedImpact);
            }
            snv.mixedGermlineImpact(mnv.mixedGermlineImpact());
        }
    }

    private static boolean isPassingSnv(SageVariant variant) {
        return variant.isPassing() && variant.isSnv();
    }

    private static boolean isMixedGermlineMnv(SageVariant variant) {
        return variant.filters().contains(SoftFilter.MIXED_GERMLINE_SOMATIC_MNV.toString()) && variant.isMnv();
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
