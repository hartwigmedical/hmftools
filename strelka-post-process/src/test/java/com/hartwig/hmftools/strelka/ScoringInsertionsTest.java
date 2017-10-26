package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.TestUtils.buildSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.scores.ImmutableReadScore;
import com.hartwig.hmftools.strelka.scores.ReadScore;
import com.hartwig.hmftools.strelka.scores.ReadType;

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
        final SAMRecord reference = buildSamRecord(1, "9M", "GATCCGATC", false);
        final ReadScore score = Scoring.getReadScore(reference, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumor() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I7M", "GATCTCCGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSandDEL() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I1M2D4M", "GATCTCCGA", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithDELPre() {
        final SAMRecord tumor = buildSamRecord(1, "1D1M2I7M", "ATCTCCGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithDELAfter() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I1M2D4M", "GATCTGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSPre() {
        final SAMRecord tumor = buildSamRecord(1, "1M3I1M2I7M", "GGCCATCTCCGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSinTumorWithINSAfter() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I2M4I5M", "GATCTCAAAACGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotdetectINSinRefWithSNV() {
        final SAMRecord tumor = buildSamRecord(1, "9M", "GAGGCGATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectINSinTumorWithShortRead() {
        final SAMRecord referenceShortRead = buildSamRecord(1, "3M", "GAT", false);
        final ReadScore score = Scoring.getReadScore(referenceShortRead, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void detectsINSatEndOfRead() {
        final SAMRecord tumor = buildSamRecord(1, "2M2I", "GATC", false);
        final ReadScore score = Scoring.getReadScore(tumor, INSERTION);
        assertEquals(ReadType.ALT, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void doesNotDetectINSatEndOfReadInRef() {
        final SAMRecord reference = buildSamRecord(1, "4M", "GATC", false);
        final ReadScore score = Scoring.getReadScore(reference, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertTrue(score.score() > 0);
    }

    @Test
    public void computesScoreForINSinRef() {
        final SAMRecord ref = buildSamRecord(2, "2M", "AA", "++", false);
        assertEquals(ImmutableReadScore.of(ReadType.REF, 10), Scoring.getReadScore(ref, INSERTION));
    }

    @Test
    public void computesScoreForINSinTumor() {
        //MIVO: insertion with qualities 30, 20, 40 --> average = 30
        final SAMRecord alt = buildSamRecord(2, "1M2I", "ATC", "?5I", false);
        assertEquals(ImmutableReadScore.of(ReadType.ALT, 30), Scoring.getReadScore(alt, INSERTION));
    }

    @Test
    public void computesScoreForINSinOther() {
        final SAMRecord otherSNV = buildSamRecord(2, "2M", "CC", "C", false);
        assertEquals(ImmutableReadScore.of(ReadType.REF, 34), Scoring.getReadScore(otherSNV, INSERTION));
    }

    @Test
    public void computesScoreForINSinReadWithDeletionOnVariantPos() {
        final SAMRecord deleted = buildSamRecord(1, "1M2D1M", "AA", "FD", false);
        assertEquals(ImmutableReadScore.of(ReadType.MISSING, 0), Scoring.getReadScore(deleted, INSERTION));
    }

    @Test
    public void doesNotComputeScoreForShorterINSinTumor() {
        final SAMRecord alt = buildSamRecord(2, "1M1I1M", "ATC", "?5I", false);
        assertEquals(ImmutableReadScore.of(ReadType.REF, 30), Scoring.getReadScore(alt, INSERTION));
    }

    @Test
    public void doesNotComputeScoreForLongerINSinTumor() {
        final SAMRecord alt = buildSamRecord(2, "1M3I", "ATCC", "?5II", false);
        final ReadScore score = Scoring.getReadScore(alt, INSERTION);
        assertEquals(ReadType.REF, score.type());
        assertEquals(30, score.score());
    }

    @Test
    public void doesNotComputeScoreForLongerINSandMatchInTumor() {
        final SAMRecord alt = buildSamRecord(2, "1M3I1M", "ATCCM", "?5II?", false);
        final ReadScore score = Scoring.getReadScore(alt, INSERTION);
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
        assertEquals(32, Scoring.getReadScore(cleanRecord, INSERTION).score());
        assertEquals(ReadType.REF, Scoring.getReadScore(cleanRecord, INSERTION).type());
        assertEquals(32, Scoring.getReadScore(cleanRecord, SNV).score());
        assertEquals(ReadType.REF, Scoring.getReadScore(cleanRecord, INSERTION).type());
        assertEquals(32, Scoring.getReadScore(insertMNVRecord, INSERTION).score());
        assertEquals(ReadType.ALT, Scoring.getReadScore(insertMNVRecord, INSERTION).type());
        assertEquals(32, Scoring.getReadScore(insertMNVRecord, SNV).score());
        assertEquals(ReadType.ALT, Scoring.getReadScore(insertMNVRecord, INSERTION).type());
    }

    @NotNull
    private List<SAMRecord> insertionMNVRecords() {
        final SAMRecord insertionMNVRecord = buildSamRecord(170755676, "100M1I50M",
                "AGGTTTTAGAAACAGGCTGTAAACCAGAGGGGAATAACCCACTGTGGCTAAAAAAGTAAACATAAACTTTGCAGACTCATTGGAAACAGTGTGTGGATATAAAAAAAAATATGTACTGTAACCAGAGAAAAAGAAGGCTACTTAAGAGAGA",
                false);
        final SAMRecord noMNVRecord = buildSamRecord(170755701, "151M",
                "AGTGGGGAATAACCCACTGTGGCTAAAAAATTAAACATAAACTTTTCAGACTCATTGGAAACAGTGTGTGGATATGAAAAAAATATGTACTGTAACCAGAGAAAAAGAAGGCTACTTAAGAGATAACTCTCCAATGAGTATCTTTTGTTTC",
                false);
        return Lists.newArrayList(noMNVRecord, insertionMNVRecord);
    }
}