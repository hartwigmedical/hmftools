package com.hartwig.hmftools.sage.phase;

import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.select.TierSelector;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;

class SnvSnvMerge implements Consumer<SageVariant> {

    private static final int BUFFER = 2;

    private final boolean enabled;
    private final MnvFactory factory;
    private final Consumer<SageVariant> consumer;
    private final List<SageVariant> list = Lists.newLinkedList();
    private final TierSelector tierSelector;

    SnvSnvMerge(final SageConfig config, @NotNull final Consumer<SageVariant> consumer, @NotNull final MnvFactory factory,
            @NotNull final List<GenomeRegion> panel, @NotNull final List<VariantHotspot> hotspots) {
        this.enabled = config.mnvDetection();
        this.consumer = consumer;
        this.factory = factory;
        this.tierSelector = new TierSelector(panel, hotspots);
    }

    @Override
    public void accept(@NotNull final SageVariant newEntry) {
        flush(newEntry);
        if (enabled && isPhasedSnv(newEntry)) {

            for (int i = 0; i < list.size(); i++) {
                final SageVariant oldEntry = list.get(i);
                if (isMnv(oldEntry, newEntry)) {
                    final VariantHotspot candidate = factory.merge(oldEntry.primaryTumor(), newEntry.primaryTumor());

                    boolean candidateIsHotspot = tierSelector.isHotspot(candidate);
                    boolean bothEntriesPass = newEntry.isPassing() && oldEntry.isPassing();
                    boolean onePassingOneGermlineFiltered = onePassingOneGermlineFiltered(oldEntry, newEntry);

                    if (bothEntriesPass || candidateIsHotspot || onePassingOneGermlineFiltered) {
                        SageVariant mnv = factory.mnv(newEntry.localPhaseSet(), candidate);
                        if (isPassingWithNoSupportInNormal(mnv) ||  candidateIsHotspot) {
                            newEntry.filters().add(SageVCF.MERGE_FILTER);
                            oldEntry.filters().add(SageVCF.MERGE_FILTER);
                            if (oldEntry.isSynthetic()) {
                                list.set(i, mnv);
                            } else {
                                list.add(i, mnv);
                                i++;
                            }
                        } else if (onePassingOneGermlineFiltered) {
                            mnv.filters().add(SageVCF.GERMLINE_MVN);
                            list.add(i, mnv);
                            i++;
                        }
                    }
                }
            }
        }

        list.add(newEntry);
    }

    private boolean onePassingOneGermlineFiltered(SageVariant oldVariant, SageVariant newVariant) {
        return oldVariant.isPassing() && SoftFilter.isGermlineAndNotTumorFiltered(newVariant.filters())
                || newVariant.isPassing() && SoftFilter.isGermlineAndNotTumorFiltered(oldVariant.filters());
    }

    private boolean isPassingWithNoSupportInNormal(@NotNull final SageVariant mnv) {
        return mnv.isPassing() && mnv.normal().primaryReadContext().altSupport() == 0;
    }

    private void flush(@NotNull final GenomePosition position) {
        final Iterator<SageVariant> iterator = list.iterator();
        while (iterator.hasNext()) {
            final SageVariant entry = iterator.next();
            long entryEnd = entry.position() + entry.normal().ref().length() - 1;
            if (!entry.chromosome().equals(position.chromosome()) || entryEnd < position.position() - BUFFER) {
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

    private boolean isPhasedSnv(@NotNull final SageVariant newEntry) {
        return newEntry.localPhaseSet() > 0 && !newEntry.isIndel() && !newEntry.filters().contains(SageVCF.MERGE_FILTER);
    }

    private boolean isMnv(@NotNull final SageVariant existingEntry, @NotNull final SageVariant newEntry) {
        long existingEntryEndPosition = existingEntry.position() + existingEntry.normal().ref().length() - 1;
        return isPhasedSnv(existingEntry) && existingEntry.localPhaseSet() == newEntry.localPhaseSet()
                && newEntry.position() - existingEntryEndPosition <= BUFFER;
    }
}
