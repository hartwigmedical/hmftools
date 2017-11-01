package com.hartwig.hmftools.strelka.mnv.scores;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.mnv.TestUtils;
import com.hartwig.hmftools.strelka.mnv.scores.ImmutableVariantScore;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ScoringInsertionsTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());
    private static final VariantContext INSERTION = VARIANTS.get(0);

    @Test
    public void doesNotDetectINSinRef() {
        final SAMRecord reference = TestUtils.buildSamRecord(1, "9M", "GATCCGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(reference, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumor() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I7M", "GATCTCCGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSandDEL() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I1M2D4M", "GATCTCCGA");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithDELPre() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "1D1M2I7M", "ATCTCCGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithDELAfter() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I1M2D4M", "GATCTGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSPre() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "1M3I1M2I7M", "GGCCATCTCCGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSAfter() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I2M4I5M", "GATCTCAAAACGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotdetectINSinRefWithSNV() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "9M", "GAGGCGATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectINSinTumorWithShortRead() {
        final SAMRecord referenceShortRead = TestUtils.buildSamRecord(1, "3M", "GAT");
        final VariantScore score = SamRecordScoring.getVariantScore(referenceShortRead, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSatEndOfRead() {
        final SAMRecord tumor = TestUtils.buildSamRecord(1, "2M2I", "GATC");
        final VariantScore score = SamRecordScoring.getVariantScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectINSatEndOfReadInRef() {
        final SAMRecord reference = TestUtils.buildSamRecord(1, "4M", "GATC");
        final VariantScore score = SamRecordScoring.getVariantScore(reference, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void computesScoreForINSinRef() {
        final SAMRecord ref = TestUtils.buildSamRecord(2, "2M", "AA", "++", false);
        assertEquals(ImmutableVariantScore.of(ReadType.REF, 10), SamRecordScoring.getVariantScore(ref, INSERTION));
    }

    @Test
    public void computesScoreForINSinTumor() {
        //MIVO: insertion with qualities 30, 20, 40 --> average = 30
        final SAMRecord alt = TestUtils.buildSamRecord(2, "1M2I", "ATC", "?5I", false);
        assertEquals(ImmutableVariantScore.of(ReadType.ALT, 30), SamRecordScoring.getVariantScore(alt, INSERTION));
    }

    @Test
    public void computesScoreForINSinOther() {
        final SAMRecord otherSNV = TestUtils.buildSamRecord(2, "2M", "CC", "C", false);
        assertEquals(ImmutableVariantScore.of(ReadType.REF, 34), SamRecordScoring.getVariantScore(otherSNV, INSERTION));
    }

    @Test
    public void computesScoreForINSinReadWithDeletionOnVariantPos() {
        final SAMRecord deleted = TestUtils.buildSamRecord(1, "1M2D1M", "AA", "FD", false);
        assertEquals(ImmutableVariantScore.of(ReadType.MISSING, 0), SamRecordScoring.getVariantScore(deleted, INSERTION));
    }

    @Test
    public void doesNotComputeScoreForShorterINSinTumor() {
        final SAMRecord alt = TestUtils.buildSamRecord(2, "1M1I1M", "ATC", "?5I", false);
        assertEquals(ImmutableVariantScore.of(ReadType.REF, 30), SamRecordScoring.getVariantScore(alt, INSERTION));
    }

    @Test
    public void doesNotComputeScoreForLongerINSinTumor() {
        final SAMRecord alt = TestUtils.buildSamRecord(2, "1M3I", "ATCC", "?5II", false);
        final VariantScore score = SamRecordScoring.getVariantScore(alt, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertEquals(30, score.score());
    }

    @Test
    public void doesNotComputeScoreForLongerINSandMatchInTumor() {
        final SAMRecord alt = TestUtils.buildSamRecord(2, "1M3I1M", "ATCCM", "?5II?", false);
        final VariantScore score = SamRecordScoring.getVariantScore(alt, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertEquals(30, score.score());
    }

    @Test
    public void containsInsertionMNVTest() {
        final VariantContext INSERTION = VARIANTS.get(5);
        final VariantContext SNV = VARIANTS.get(6);
        final List<SAMRecord> records = insertionMNVRecords();
        final SAMRecord cleanRecord = records.get(0);
        final SAMRecord insertMNVRecord = records.get(1);
        assertEquals(32, SamRecordScoring.getVariantScore(cleanRecord, INSERTION).score());
        assertEquals(ReadType.REF, SamRecordScoring.getVariantScore(cleanRecord, INSERTION).type());
        assertEquals(32, SamRecordScoring.getVariantScore(cleanRecord, SNV).score());
        assertEquals(ReadType.REF, SamRecordScoring.getVariantScore(cleanRecord, INSERTION).type());
        assertEquals(32, SamRecordScoring.getVariantScore(insertMNVRecord, INSERTION).score());
        assertEquals(ReadType.ALT, SamRecordScoring.getVariantScore(insertMNVRecord, INSERTION).type());
        assertEquals(32, SamRecordScoring.getVariantScore(insertMNVRecord, SNV).score());
        assertEquals(ReadType.ALT, SamRecordScoring.getVariantScore(insertMNVRecord, INSERTION).type());
    }

    @NotNull
    private List<SAMRecord> insertionMNVRecords() {
        final SAMRecord insertionMNVRecord = TestUtils.buildSamRecord(170755676, "100M1I50M",
                "AGGTTTTAGAAACAGGCTGTAAACCAGAGGGGAATAACCCACTGTGGCTAAAAAAGTAAACATAAACTTTGCAGACTCATTGGAAACAGTGTGTGGATATAAAAAAAAATATGTACTGTAACCAGAGAAAAAGAAGGCTACTTAAGAGAGA");
        final SAMRecord noMNVRecord = TestUtils.buildSamRecord(170755701, "151M",
                "AGTGGGGAATAACCCACTGTGGCTAAAAAATTAAACATAAACTTTTCAGACTCATTGGAAACAGTGTGTGGATATGAAAAAAATATGTACTGTAACCAGAGAAAAAGAAGGCTACTTAAGAGATAACTCTCCAATGAGTATCTTTTGTTTC");
        return Lists.newArrayList(noMNVRecord, insertionMNVRecord);
    }
}