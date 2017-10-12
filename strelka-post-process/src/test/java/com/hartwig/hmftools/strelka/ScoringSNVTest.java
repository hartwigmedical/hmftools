package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.MNVDetectorTest.buildSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.scores.ImmutableReadScore;
import com.hartwig.hmftools.strelka.scores.ReadScore;
import com.hartwig.hmftools.strelka.scores.ReadType;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ScoringSNVTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());
    private static final VariantContext SNV = VARIANTS.get(2);

    @Test
    public void doesNotDetectSNVinRef() {
        final SAMRecord reference = buildSamRecord(1, "6M", "GATCCG", false);
        final ReadScore score = Scoring.getReadScore(reference, SNV);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectSNVWithDeletionOnPos() {
        final SAMRecord recordWithDeletion = buildSamRecord(1, "3M1D2M", "GATCG", false);
        final ReadScore score = Scoring.getReadScore(recordWithDeletion, SNV);
        assertEquals(ReadType.MISSING, score.type());
        assertEquals(0, score.score());
    }

    @Test
    public void detectsSNVinTumor() {
        final SAMRecord tumor = buildSamRecord(1, "6M", "GATTCG", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithDELPre() {
        final SAMRecord tumor = buildSamRecord(1, "2M1D3M", "GATCG", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithDELAfter() {
        final SAMRecord tumor = buildSamRecord(1, "4M1D1M", "GATTG", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithINSPre() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I4M", "GACCTTCG", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithINSAfter() {
        final SAMRecord tumor = buildSamRecord(1, "4M2I2M", "GATTAACG", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVatEndOfRead() {
        final SAMRecord tumor = buildSamRecord(1, "4M", "GATT", false);
        final ReadScore score = Scoring.getReadScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectSNVatEndOfReadInRef() {
        final SAMRecord reference = buildSamRecord(1, "4M", "GATC", false);
        final ReadScore score = Scoring.getReadScore(reference, SNV);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void computesScoreForSNVinRef() {
        final SAMRecord ref = buildSamRecord(4, "1M", "C", "+", false);
        assertEquals(ImmutableReadScore.of(ReadType.REF, 10), Scoring.getReadScore(ref, SNV));
    }

    @Test
    public void computesScoreForSNVinTumor() {
        final SAMRecord alt = buildSamRecord(4, "1M", "T", "P", false);
        assertEquals(ImmutableReadScore.of(ReadType.ALT, 47), Scoring.getReadScore(alt, SNV));
    }

    @Test
    public void computesScoreForSNVinOther() {
        final SAMRecord otherSNV = buildSamRecord(4, "1M", "A", "F", false);
        assertEquals(ImmutableReadScore.of(ReadType.REF, 37), Scoring.getReadScore(otherSNV, SNV));
    }

    @Test
    public void computesScoreForSNVinReadWithDeletionOnVariantPos() {
        final SAMRecord deleted = buildSamRecord(3, "1M1D1M", "AA", "FD", false);
        assertEquals(ImmutableReadScore.of(ReadType.MISSING, 0), Scoring.getReadScore(deleted, SNV));
    }
}
