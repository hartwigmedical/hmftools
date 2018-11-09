package com.hartwig.hmftools.common.hotspot;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class HotspotEvidenceFactory {

    private final LinkedHashSet<VariantHotspot> hotspots;

    public HotspotEvidenceFactory(final ListMultimap<Chromosome, VariantHotspot> hotspots) {
        this.hotspots = new LinkedHashSet<>(hotspots.values());
    }

    @NotNull
    public List<HotspotEvidence> evidence(@NotNull final List<Pileup> tumor, @NotNull final List<Pileup> normal) {
        final List<HotspotEvidence> result = Lists.newArrayList();

        result.addAll(hotspotEvidence(tumor, normal));
        result.addAll(inframeIndelEvidence(tumor, normal));

        Collections.sort(result);
        return result;
    }

    @NotNull
    private List<HotspotEvidence> hotspotEvidence(@NotNull final List<Pileup> tumor, @NotNull final List<Pileup> normal) {
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
                    result.add(fromHotspot(hotspot, tumorPileup, optionalNormalPileup));
                }
            }
        }

        return result;
    }

    @NotNull
    static HotspotEvidence fromHotspot(@NotNull final VariantHotspot hotspot, @NotNull final Pileup tumor, @NotNull final Optional<Pileup> normal) {
        int tumorEvidence = evidence(tumor, hotspot);
        assert (tumorEvidence > 0);
        return ImmutableHotspotEvidence.builder()
                .from(hotspot)
                .type(HotspotEvidenceType.fromVariantHotspot(hotspot))
                .alt(hotspot.alt())
                .ref(hotspot.ref())
                .qualityScore(qualityScore(tumor, hotspot))
                .tumorEvidence(tumorEvidence)
                .tumorReads(tumor.readCount())
                .normalEvidence(normal.map(x -> evidence(x, hotspot)).orElse(0))
                .normalReads(normal.map(Pileup::readCount).orElse(0))
                .build();
    }

    @NotNull
    private List<HotspotEvidence> inframeIndelEvidence(@NotNull final List<Pileup> tumor, @NotNull final List<Pileup> normal) {
        final List<HotspotEvidence> result = Lists.newArrayList();

        final GenomePositionSelector<Pileup> normalSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(normal));
        for (Pileup pileup : tumor) {
            result.addAll(inframeIndelEvidence(pileup,
                    normalSelector.select(GenomePositions.create(pileup.chromosome(), pileup.position()))));

        }

        return result;
    }

    @NotNull
    static List<HotspotEvidence> inframeIndelEvidence(@NotNull final Pileup tumor, @NotNull final Optional<Pileup> normal) {
        final List<HotspotEvidence> result = Lists.newArrayList();

        final ImmutableHotspotEvidence.Builder builder = ImmutableHotspotEvidence.builder()
                .chromosome(tumor.chromosome())
                .type(HotspotEvidenceType.INFRAME)
                .position(tumor.position());

        for (final String insert : tumor.inframeInserts()) {
            result.add(builder.ref(tumor.referenceBase())
                    .alt(insert)
                    .qualityScore(tumor.insertScore(insert))
                    .tumorEvidence(tumor.insertCount(insert))
                    .tumorReads(tumor.readCount())
                    .normalEvidence(normal.map(x -> x.insertCount(insert)).orElse(0))
                    .normalReads(normal.map(Pileup::readCount).orElse(0))
                    .build());
        }

        for (final String del : tumor.inframeDeletes()) {
            result.add(builder.alt(tumor.referenceBase())
                    .ref(del)
                    .qualityScore(tumor.deleteScore(del))
                    .tumorEvidence(tumor.deleteCount(del))
                    .tumorReads(tumor.readCount())
                    .normalEvidence(normal.map(x -> x.deleteCount(del)).orElse(0))
                    .normalReads(normal.map(Pileup::readCount).orElse(0))
                    .build());
        }

        return result;
    }

    private static int evidence(@NotNull final Pileup pileup, @NotNull final VariantHotspot hotspot) {

        if (hotspot.isSNV()) {
            return pileup.mismatchCount(hotspot.alt().charAt(0));
        }

        if (hotspot.isInsert()) {
            return pileup.insertCount(hotspot.alt());
        }

        if (hotspot.isDelete()) {
            return pileup.deleteCount(hotspot.ref());
        }

        return 0;
    }

    @VisibleForTesting
    static int qualityScore(@NotNull final Pileup tumor, @NotNull final VariantHotspot hotspot) {
        if (hotspot.isSNV()) {
            return tumor.mismatchScore(hotspot.alt().charAt(0));
        }

        if (hotspot.isInsert()) {
            return tumor.insertScore(hotspot.alt());
        }

        if (hotspot.isDelete()) {
            return tumor.deleteScore(hotspot.ref());
        }

        return 0;
    }
}
