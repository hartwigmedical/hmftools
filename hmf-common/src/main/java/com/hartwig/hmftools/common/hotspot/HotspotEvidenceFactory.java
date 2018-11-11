package com.hartwig.hmftools.common.hotspot;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class HotspotEvidenceFactory {

    private final Collection<VariantHotspot> hotspots;

    public HotspotEvidenceFactory(final ListMultimap<Chromosome, VariantHotspot> hotspots) {
        this.hotspots = hotspots.values();
    }

    @NotNull
    public List<HotspotEvidence> evidence(@NotNull final List<Pileup> tumor, @NotNull final List<Pileup> normal) {

        final Map<VariantHotspot, HotspotEvidence> resultMap = Maps.newHashMap();
        for (HotspotEvidence evidence : hotspotEvidence(tumor, normal)) {
            resultMap.put(fromEvidence(evidence), evidence);
        }

        for (HotspotEvidence evidence : inframeIndelEvidence(tumor, normal)) {
            resultMap.putIfAbsent(fromEvidence(evidence), evidence);
        }

        final List<HotspotEvidence> result = Lists.newArrayList(resultMap.values());
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

                if (hotspot.isMNV()) {
                    final GenomeRegion mnvRegion = GenomeRegionFactory.create(hotspot.chromosome(),
                            hotspot.position(),
                            hotspot.position() + hotspot.ref().length() - 1);

                    final MNVEvidence tumorMnvEvidence = new MNVEvidence(hotspot);
                    tumorSelector.select(mnvRegion, tumorMnvEvidence);
                    if (tumorMnvEvidence.evidence() > 0) {
                        final MNVEvidence normalMnvEvidence = new MNVEvidence(hotspot);
                        tumorSelector.select(mnvRegion, normalMnvEvidence);
                        result.add(fromMNV(hotspot, tumorMnvEvidence, normalMnvEvidence));
                    }
                } else {
                    final Pileup tumorPileup = optionalTumorPileup.get();
                    int tumorEvidence = evidence(tumorPileup, hotspot);
                    if (tumorEvidence > 0) {
                        final Optional<Pileup> optionalNormalPileup = normalSelector.select(hotspot);
                        result.add(fromHotspot(tumorEvidence, hotspot, tumorPileup, optionalNormalPileup));
                    }
                }
            }
        }

        return result;
    }

    @NotNull
    private static HotspotEvidence fromHotspot(int tumorEvidence, @NotNull final VariantHotspot hotspot, @NotNull final Pileup tumor,
            @NotNull final Optional<Pileup> normal) {
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
    private static HotspotEvidence fromMNV(@NotNull final VariantHotspot hotspot, @NotNull final MNVEvidence tumor,
            @NotNull final MNVEvidence normal) {
        return ImmutableHotspotEvidence.builder()
                .from(hotspot)
                .type(HotspotEvidenceType.MNV)
                .alt(hotspot.alt())
                .ref(hotspot.ref())
                .qualityScore(tumor.score())
                .tumorEvidence(tumor.evidence())
                .tumorReads(tumor.reads())
                .normalEvidence(normal.evidence())
                .normalReads(normal.reads())
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
    private static List<HotspotEvidence> inframeIndelEvidence(@NotNull final Pileup tumor, @NotNull final Optional<Pileup> normal) {
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
    private static int qualityScore(@NotNull final Pileup tumor, @NotNull final VariantHotspot hotspot) {
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

    @NotNull
    private static VariantHotspot fromEvidence(@NotNull final HotspotEvidence evidence) {
        return ImmutableVariantHotspot.builder()
                .chromosome(evidence.chromosome())
                .position(evidence.position())
                .ref(evidence.ref())
                .alt(evidence.alt())
                .build();
    }
}
