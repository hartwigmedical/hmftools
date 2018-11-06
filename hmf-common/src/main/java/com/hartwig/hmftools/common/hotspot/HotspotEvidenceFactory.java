package com.hartwig.hmftools.common.hotspot;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public class HotspotEvidenceFactory {

    private final LinkedHashSet<VariantHotspot> hotspots;

    public HotspotEvidenceFactory(final ListMultimap<Chromosome, VariantHotspot> hotspots) {
        this.hotspots = new LinkedHashSet<>(hotspots.values());
    }

    @NotNull
    public List<HotspotEvidence> hotspotEvidence(@NotNull final List<Pileup> tumor, @NotNull final List<Pileup> normal) {
        final List<HotspotEvidence> result = Lists.newArrayList();

        final GenomePositionSelector<Pileup> tumorSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(tumor));
        final GenomePositionSelector<Pileup> normalSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(normal));

        for (final VariantHotspot hotspot : hotspots) {

            final Optional<Pileup> optionalTumorPileup = tumorSelector.select(hotspot);
            if (optionalTumorPileup.isPresent()) {

                final Pileup tumorPileup = optionalTumorPileup.get();
                int tumorEvidence = evidence(tumorPileup, hotspot);
                if (tumorEvidence > 0) {
                    final Optional<Pileup> optionalNormalPileup = normalSelector.select(hotspot);
                    result.add(ImmutableHotspotEvidence.builder()
                            .from(hotspot)
                            .type(HotspotEvidenceType.fromVariantHotspot(hotspot))
                            .alt(hotspot.alt())
                            .ref(hotspot.ref())
                            .tumorEvidence(tumorEvidence)
                            .tumorReads(tumorPileup.readCount())
                            .normalEvidence(optionalNormalPileup.map(x -> evidence(x, hotspot)).orElse(0))
                            .normalReads(optionalNormalPileup.map(Pileup::readCount).orElse(0))
                            .build());
                }
            }
        }
        return result;
    }

    private int evidence(@NotNull final Pileup pileup, @NotNull VariantHotspot hotspot) {

        if (hotspot.isSNV()) {
            return pileup.mismatchCount(hotspot.alt().charAt(0));
        }

        if (hotspot.isInsert()) {
            return pileup.insertionCounts().getOrDefault(hotspot.alt(), 0);
        }

        if (hotspot.isDelete()) {
            return pileup.deletionCounts().getOrDefault(hotspot.ref(), 0);
        }

        return 0;
    }
}
