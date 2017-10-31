package com.hartwig.hmftools.strelka;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.strelka.scores.ReadType;
import com.hartwig.hmftools.strelka.scores.VariantScore;
import com.hartwig.hmftools.strelka.scores.VariantScoreAggregate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MNVScore {
    private static final double MNV_SCORE_THRESHOLD = 0.8;

    @NotNull
    public abstract Map<VariantContext, VariantScoreAggregate> aggregatedScores();

    @NotNull
    public abstract Map<VariantContext, Integer> missingPerPosition();

    public abstract long scoreWithOneAlt();

    public abstract long scoreWithAllAlts();

    //    @NotNull
    //    @Value.Derived
    //    public List<VariantContext> variants() {
    //        return aggregatedScores().keySet().stream().sorted(Comparator.comparing(VariantContext::getStart)).collect(Collectors.toList());
    //    }

    @NotNull
    static MNVScore addReadScore(@NotNull final MNVScore mnvScore, @NotNull final Map<VariantContext, VariantScore> variantScores) {
        final Map<VariantContext, VariantScoreAggregate> aggregatedScores = Maps.newHashMap();
        final Map<VariantContext, Integer> missingCounts = Maps.newHashMap();
        boolean oneAlt = false;
        boolean allAlts = true;
        int scoreTotal = 0;
        for (final VariantContext variant : mnvScore.aggregatedScores().keySet()) {
            final VariantScore variantScore = variantScores.get(variant);
            scoreTotal += variantScore.score();
            aggregatedScores.put(variant, mnvScore.aggregatedScores().get(variant).addScore(variantScore));
            if (variantScore.type() != ReadType.ALT) {
                allAlts = false;
            }
            if (variantScore.type() == ReadType.ALT) {
                oneAlt = true;
            }
            if (variantScore.type() == ReadType.MISSING) {
                missingCounts.put(variant, mnvScore.missingPerPosition().get(variant) + 1);
            } else {
                missingCounts.put(variant, mnvScore.missingPerPosition().get(variant));
            }
        }
        if (allAlts) {
            return ImmutableMNVScore.of(aggregatedScores, missingCounts, mnvScore.scoreWithOneAlt() + scoreTotal,
                    mnvScore.scoreWithAllAlts() + scoreTotal);
        } else if (oneAlt) {
            return ImmutableMNVScore.of(aggregatedScores, missingCounts, mnvScore.scoreWithOneAlt() + scoreTotal,
                    mnvScore.scoreWithAllAlts());
        } else {
            return ImmutableMNVScore.of(aggregatedScores, mnvScore.missingPerPosition(), mnvScore.scoreWithOneAlt(),
                    mnvScore.scoreWithAllAlts());
        }
    }

    double frequency() {
        final long correction = missingPerPosition().entrySet().stream().mapToInt(entry -> {
            final int averagePerPosition = aggregatedScores().get(entry.getKey()).average();
            return averagePerPosition * missingPerPosition().get(entry.getKey());
        }).sum();
        final long correctedScoreWithOneVariant = scoreWithOneAlt() + correction;
        if (correctedScoreWithOneVariant == 0) {
            return 0;
        } else {
            return (double) scoreWithAllAlts() / correctedScoreWithOneVariant;
        }
    }

    boolean isMNV() {
        return frequency() > MNV_SCORE_THRESHOLD;
    }

    @NotNull
    public static MNVScore of(@NotNull final List<VariantContext> variants) {
        final Map<VariantContext, VariantScoreAggregate> scores =
                variants.stream().collect(Collectors.toMap(Function.identity(), (variant) -> VariantScoreAggregate.newScore()));
        final Map<VariantContext, Integer> missing = variants.stream().collect(Collectors.toMap(Function.identity(), (variant) -> 0));
        return ImmutableMNVScore.of(scores, missing, 0L, 0L);
    }
}
