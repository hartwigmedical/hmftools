package com.hartwig.hmftools.strelka;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.strelka.scores.ReadScore;
import com.hartwig.hmftools.strelka.scores.ReadType;
import com.hartwig.hmftools.strelka.scores.VariantScore;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MNVScore {
    private static double MNV_SCORE_THRESHOLD = 0.8;

    @NotNull
    public abstract Map<VariantContext, VariantScore> scoresPerPosition();

    @NotNull
    public abstract Map<VariantContext, Integer> missingPerPosition();

    public abstract long scoreWithOneAlt();

    public abstract long scoreWithAllAlts();

    static MNVScore addReadScore(@NotNull final MNVScore mnvScore, @NotNull final Map<VariantContext, ReadScore> readScore) {
        final Map<VariantContext, VariantScore> updatedScores = Maps.newHashMap();
        final Map<VariantContext, Integer> missingCounts = Maps.newHashMap();
        boolean oneAlt = false;
        boolean allAlts = true;
        int scoreSum = 0;
        for (final VariantContext variant : mnvScore.scoresPerPosition().keySet()) {
            final ReadScore readScoreAtVariant = readScore.get(variant);
            scoreSum += readScoreAtVariant.score();
            updatedScores.put(variant, mnvScore.scoresPerPosition().get(variant).addScore(readScoreAtVariant));
            if (readScoreAtVariant.type() != ReadType.ALT) {
                allAlts = false;
            }
            if (readScoreAtVariant.type() == ReadType.ALT) {
                oneAlt = true;
            }
            if (readScoreAtVariant.type() == ReadType.MISSING) {
                missingCounts.put(variant, mnvScore.missingPerPosition().get(variant) + 1);
            } else {
                missingCounts.put(variant, mnvScore.missingPerPosition().get(variant));
            }
        }
        if (allAlts) {
            return ImmutableMNVScore.of(updatedScores, missingCounts, mnvScore.scoreWithOneAlt() + scoreSum,
                    mnvScore.scoreWithAllAlts() + scoreSum);
        } else if (oneAlt) {
            return ImmutableMNVScore.of(updatedScores, missingCounts, mnvScore.scoreWithOneAlt() + scoreSum, mnvScore.scoreWithAllAlts());
        } else {
            return ImmutableMNVScore.of(updatedScores, mnvScore.missingPerPosition(), mnvScore.scoreWithOneAlt(),
                    mnvScore.scoreWithAllAlts());
        }
    }

    double frequency() {
        final long correction = missingPerPosition().entrySet().stream().mapToInt(entry -> {
            final int averagePerPosition = scoresPerPosition().get(entry.getKey()).average();
            return averagePerPosition * missingPerPosition().get(entry.getKey());
        }).sum();
        final long correctedScoreWithOneVariant = scoreWithOneAlt() + correction;
        return (double) scoreWithAllAlts() / correctedScoreWithOneVariant;
    }

    boolean isMNV() {
        return frequency() > MNV_SCORE_THRESHOLD;
    }

    public static MNVScore of(@NotNull final List<VariantContext> variants) {
        final Map<VariantContext, VariantScore> scores =
                variants.stream().collect(Collectors.toMap(Function.identity(), (variant) -> VariantScore.newScore()));
        final Map<VariantContext, Integer> missing = variants.stream().collect(Collectors.toMap(Function.identity(), (variant) -> 0));
        return ImmutableMNVScore.of(scores, missing, 0L, 0L);
    }
}
