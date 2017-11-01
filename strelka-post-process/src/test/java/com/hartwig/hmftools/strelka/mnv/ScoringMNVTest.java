package com.hartwig.hmftools.strelka.mnv;

import static com.hartwig.hmftools.strelka.mnv.TestUtils.build2VariantScores;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.mnv.scores.ImmutableVariantScore;
import com.hartwig.hmftools.strelka.mnv.scores.ReadType;
import com.hartwig.hmftools.strelka.mnv.scores.VariantScore;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ScoringMNVTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    @Test
    public void correctlyComputes100PercentFrequency() {
        final List<Pair<VariantScore, VariantScore>> mnvScores = Lists.newArrayList();
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.REF, 10), ImmutableVariantScore.of(ReadType.REF, 20)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15)));
        final MNVScore scores = build2VariantScores(Lists.newArrayList(VARIANTS.get(0), VARIANTS.get(1)), mnvScores);
        assertEquals(1.0, scores.frequency(), 0.000001);
    }

    @Test
    public void correctlyComputes50PercentFrequency() {
        final List<Pair<VariantScore, VariantScore>> mnvScores = Lists.newArrayList();
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 10), ImmutableVariantScore.of(ReadType.REF, 20)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15)));
        final MNVScore scores = build2VariantScores(Lists.newArrayList(VARIANTS.get(0), VARIANTS.get(1)), mnvScores);
        assertEquals(0.5, scores.frequency(), 0.000001);
    }

    @Test
    public void correctlyComputes50PercentFrequencyForMissingRead() {
        final List<Pair<VariantScore, VariantScore>> mnvScores = Lists.newArrayList();
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.MISSING, 0)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15)));
        final MNVScore scores = build2VariantScores(Lists.newArrayList(VARIANTS.get(0), VARIANTS.get(1)), mnvScores);
        assertEquals(0.5, scores.frequency(), 0.000001);
    }

    @Test
    public void correctlyComputes33PercentFrequency() {
        final List<Pair<VariantScore, VariantScore>> mnvScores = Lists.newArrayList();
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 5), ImmutableVariantScore.of(ReadType.REF, 20)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.REF, 25), ImmutableVariantScore.of(ReadType.ALT, 10)));
        final MNVScore scores = build2VariantScores(Lists.newArrayList(VARIANTS.get(0), VARIANTS.get(1)), mnvScores);
        assertEquals(0.33, scores.frequency(), 0.01);
    }

    @Test
    public void correctlyComputes80PercentFrequency() {
        final List<Pair<VariantScore, VariantScore>> mnvScores = Lists.newArrayList();
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 5), ImmutableVariantScore.of(ReadType.REF, 25)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 10), ImmutableVariantScore.of(ReadType.ALT, 20)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 13), ImmutableVariantScore.of(ReadType.ALT, 17)));
        mnvScores.add(ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 27), ImmutableVariantScore.of(ReadType.ALT, 3)));
        final MNVScore scores = build2VariantScores(Lists.newArrayList(VARIANTS.get(0), VARIANTS.get(1)), mnvScores);
        assertEquals(0.8, scores.frequency(), 0.000001);
    }
}
