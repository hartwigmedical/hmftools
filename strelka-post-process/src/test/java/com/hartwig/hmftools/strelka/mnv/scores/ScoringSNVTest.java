package com.hartwig.hmftools.strelka.mnv.scores;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.mnv.TestUtils;
import com.hartwig.hmftools.strelka.mnv.scores.ImmutableVariantScore;

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
        final SAMRecord reference = TestUtils.buildSamRecord(1, "6M", "GATCCG");
        final VariantScore score = SamRecordScoring.getVariantScore(reference, SNV);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectSNVWithDeletionOnPos() {
        final SAMRecord recordWithDeletion = TestUtils.buildSamRecord(1, "3M1D2M", "GATCG");
        final VariantScore score = SamRecordScoring.getVariantScore(recordWithDeletion, SNV);
        assertEquals(ReadType.MISSING, score.type());
        assertEquals(0, score.score());
    }

    @Test
    public void detectsSNVinTumor() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "6M", "GATTCG");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithDELPre() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M1D3M", "GATCG");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithDELAfter() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "4M1D1M", "GATTG");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithINSPre() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I4M", "GACCTTCG");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVinTumorWithINSAfter() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "4M2I2M", "GATTAACG");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsSNVatEndOfRead() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "4M", "GATT");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, SNV);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectSNVatEndOfReadInRef() {
        final SAMRecord reference = TestUtils.buildSamRecord(1, "4M", "GATC");
        final VariantScore score = SamRecordScoring.getVariantScore(reference, SNV);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void computesScoreForSNVinRef() {
        final SAMRecord ref = TestUtils.buildSamRecord(4, "1M", "C", "+", false);
        assertEquals(ImmutableVariantScore.of(ReadType.REF, 10), SamRecordScoring.getVariantScore(ref, SNV));
    }

    @Test
    public void computesScoreForSNVinTumor() {
        final SAMRecord alt = TestUtils.buildSamRecord(4, "1M", "T", "P", false);
        assertEquals(ImmutableVariantScore.of(ReadType.ALT, 47), SamRecordScoring.getVariantScore(alt, SNV));
    }

    @Test
    public void computesScoreForSNVinOther() {
        final SAMRecord otherSNV = TestUtils.buildSamRecord(4, "1M", "A", "F", false);
        assertEquals(ImmutableVariantScore.of(ReadType.REF, 37), SamRecordScoring.getVariantScore(otherSNV, SNV));
    }

    @Test
    public void computesScoreForSNVinReadWithDeletionOnVariantPos() {
        final SAMRecord deleted = TestUtils.buildSamRecord(3, "1M1D1M", "AA", "FD", false);
        assertEquals(ImmutableVariantScore.of(ReadType.MISSING, 0), SamRecordScoring.getVariantScore(deleted, SNV));
    }
}
