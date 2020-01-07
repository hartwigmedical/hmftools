package com.hartwig.hmftools.sage.phase;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import javax.annotation.concurrent.NotThreadSafe;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.select.TierSelector;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@NotThreadSafe
public class MnvCandidates implements Consumer<AltContext> {

    private static final Logger LOGGER = LogManager.getLogger(MnvCandidates.class);

    private static final int BUFFER = 2;

    private final SageConfig config;
    private final boolean enabled;
    private final TierSelector tierSelector;
    private final Consumer<AltContext> consumer;
    private final IndexedFastaSequenceFile refGenome;
    private final List<AltContext> list = Lists.newLinkedList();
    private final List<MnvCandidate> mnvCandidates = Lists.newArrayList();


    public MnvCandidates(@NotNull  final SageConfig config, @NotNull final Consumer<AltContext> consumer, @NotNull final IndexedFastaSequenceFile refGenome,
            @NotNull  final List<GenomeRegion> panel, @NotNull final List<VariantHotspot> hotspots) {
        this.enabled = config.mnvDetection();
        this.config = config;
        this.consumer = consumer;
        this.refGenome = refGenome;
        this.tierSelector = new TierSelector(panel, hotspots);
    }

    @NotNull
    public List<MnvCandidate> mvnCandidates() {
        Collections.sort(mnvCandidates);
        return mnvCandidates;
    }

    @Override
    public void accept(@NotNull final AltContext newEntry) {
        flush(newEntry);

        if (enabled && isPhasedSnv(newEntry)) {
            final List<MnvCandidate> newCandidates = Lists.newArrayList();
            boolean isNewPassing = isPassing(newEntry);

            for (final MnvCandidate oldEntry : mnvCandidates) {
                if (isMnv(oldEntry, newEntry)) {
                    final VariantHotspot mnv = createMnv(oldEntry.mnv(), newEntry);
                    if (isNewPassing || tierSelector.isHotspot(mnv)) {
                        newCandidates.add(new MnvCandidate(mnv, newEntry.localPhaseSet()));
                    }
                }
            }

            for (final AltContext oldEntry : list) {
                boolean isOldPassing = isPassing(oldEntry);

                if (isMnv(oldEntry, newEntry)) {
                    boolean bothEntriesPass = isNewPassing && isOldPassing;
                    final VariantHotspot mnv = createMnv(oldEntry, newEntry);

                    if (bothEntriesPass || tierSelector.isHotspot(mnv)) {
                        newCandidates.add(new MnvCandidate(mnv, newEntry.localPhaseSet()));
                    }
                }
            }
            mnvCandidates.addAll(newCandidates);
        }

        list.add(newEntry);
    }

    private void flush(@NotNull final GenomePosition position) {
        final Iterator<AltContext> iterator = list.iterator();
        while (iterator.hasNext()) {
            final AltContext entry = iterator.next();
            long entryEnd = entry.position() + entry.ref().length() - 1;
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

    private boolean isPhasedSnv(@NotNull final AltContext newEntry) {
        return newEntry.localPhaseSet() > 0 && newEntry.isSNV();
    }

    private boolean isMnv(@NotNull final AltContext existingEntry, @NotNull final AltContext newEntry) {
        long existingEntryEndPosition = existingEntry.position() + existingEntry.ref().length() - 1;
        return isPhasedSnv(existingEntry) && existingEntry.localPhaseSet() == newEntry.localPhaseSet()
                && newEntry.position() - existingEntryEndPosition <= BUFFER;
    }

    private boolean isMnv(@NotNull final MnvCandidate existingEntry, @NotNull final AltContext newEntry) {
        long existingEntryEndPosition = existingEntry.mnv().position() + existingEntry.mnv().ref().length() - 1;
        return existingEntry.lps() == newEntry.localPhaseSet()
                && newEntry.position() - existingEntryEndPosition <= BUFFER;
    }

    private boolean isPassing(@NotNull final AltContext newEntry) {
        final SageVariantTier tier = tierSelector.tier(newEntry);
        return config.filter().tumorFilters(tier, newEntry).isEmpty();
    }

    @NotNull
    private VariantHotspot createMnv(@NotNull final VariantHotspot left, @NotNull final AltContext right) {
        int mnvLength = (int) (right.position() - left.position() + 1);
        int additionalLength = mnvLength - left.alt().length();

        try {
            final String alt = left.alt() + right.primaryReadContext().readContext().mnvAdditionalAlt(additionalLength);
            final String ref = refGenome.getSubsequenceAt(left.chromosome(), left.position(), right.position()).getBaseString();
            return ImmutableVariantHotspotImpl.builder().from(left).ref(ref).alt(alt).build();
        } catch (Exception e) {
            LOGGER.error("Unable to merge {}:{} with {}", left.chromosome(), left.position(), right.position());
            throw e;
        }
    }

}
