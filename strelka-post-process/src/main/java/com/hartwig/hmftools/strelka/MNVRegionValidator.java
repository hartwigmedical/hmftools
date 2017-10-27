package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.SamRecordScoring.scoresPerVariant;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.strelka.scores.VariantScore;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MNVRegionValidator {
    @NotNull
    public abstract PotentialMNVRegion region();

    @NotNull
    public abstract Map<Integer, GapReads> gapReads();

    @NotNull
    public abstract Map<PotentialMNV, MNVScore> mnvScores();

    @Value.Lazy
    @NotNull
    public Map<Integer, Character> mostFrequentReads() {
        return gapReads().entrySet()
                .stream()
                .map(entry -> ImmutablePair.of(entry.getKey(), entry.getValue().mostFrequentRead()))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
    }

    @Value.Lazy
    @NotNull
    public Set<PotentialMNV> validMnvs() {
        final Set<PotentialMNV> result = Sets.newHashSet();
        final List<PotentialMNV> mnvsSortedByLength = mnvScores().entrySet()
                .stream()
                .filter(entry -> entry.getValue().isMNV())
                .sorted(mnvEntryComparator())
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
        mnvsSortedByLength.forEach(mnv -> {
            if (mnv.variants()
                    .stream()
                    .noneMatch(variant -> result.stream().anyMatch(resultMnv -> resultMnv.variants().contains(variant)))) {
                result.add(mnv);
            }
        });
        return result;
    }

    @NotNull
    public static MNVRegionValidator of(@NotNull final PotentialMNVRegion region) {
        final Map<Integer, GapReads> gapReads =
                region.gapPositions().stream().collect(Collectors.toMap(Function.identity(), (position) -> GapReads.empty()));
        Map<PotentialMNV, MNVScore> scores = region.potentialMnvs()
                .stream()
                .collect(Collectors.toMap(Function.identity(), potentialMNV -> MNVScore.of(potentialMNV.variants())));
        return ImmutableMNVRegionValidator.of(region, gapReads, scores);
    }

    @NotNull
    MNVRegionValidator addSamRecord(@NotNull final SAMRecord record) {
        if (goodRead(record) && containsAllMNVPositions(record)) {
            final Map<VariantContext, VariantScore> samRecordScores = scoresPerVariant(record, region().variantsInPotentialMnvs());
            final Map<PotentialMNV, MNVScore> scores = mnvScores().entrySet()
                    .stream()
                    .map(entry -> Pair.of(entry.getKey(), MNVScore.addReadScore(entry.getValue(), samRecordScores)))
                    .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
            final Map<Integer, GapReads> gapReads = gapReads().entrySet()
                    .stream()
                    .map(entry -> Pair.of(entry.getKey(), GapReads.addRead(entry.getValue(), getReadAtPosition(record, entry.getKey()))))
                    .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
            return ImmutableMNVRegionValidator.of(region(), gapReads, scores);
        }
        return this;
    }

    @NotNull
    private static Comparator<Map.Entry<PotentialMNV, MNVScore>> mnvEntryComparator() {
        final Comparator<Map.Entry<PotentialMNV, MNVScore>> sizeComparator =
                Comparator.comparing(mnvEntry -> mnvEntry.getKey().variants().size());
        final Comparator<Map.Entry<PotentialMNV, MNVScore>> frequencyComparator =
                Comparator.comparing(mnvEntry -> mnvEntry.getValue().frequency());
        final Comparator<Map.Entry<PotentialMNV, MNVScore>> startPositionComparator =
                Comparator.comparing(mnvEntry -> mnvEntry.getKey().start());
        return sizeComparator.reversed().thenComparing(frequencyComparator.reversed()).thenComparing(startPositionComparator);
    }

    private boolean goodRead(@NotNull final SAMRecord record) {
        return !record.getDuplicateReadFlag();
    }

    private boolean containsAllMNVPositions(@NotNull final SAMRecord record) {
        final VariantContext lastVariant = region().variants().get(region().variants().size() - 1);
        return record.getAlignmentStart() <= region().start()
                && record.getAlignmentEnd() >= lastVariant.getStart() + lastVariant.getReference().length() - 1;
    }

    @NotNull
    private static Character getReadAtPosition(@NotNull final SAMRecord record, final int position) {
        final int recordPosition = record.getReadPositionAtReferencePosition(position);
        if (recordPosition == 0) {
            return 'N';
        } else {
            return record.getReadString().charAt(recordPosition - 1);
        }
    }
}
