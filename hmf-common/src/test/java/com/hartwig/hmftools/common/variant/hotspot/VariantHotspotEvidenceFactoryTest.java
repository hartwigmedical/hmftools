package com.hartwig.hmftools.common.variant.hotspot;

import static com.hartwig.hmftools.common.utils.sam.SAMRecords.PHRED_OFFSET;
import static com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidenceFactory.create;
import static com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidenceFactory.isVariantPartOfLargerMNV;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;

import com.hartwig.hmftools.common.utils.sam.SAMRecordsTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantHotspotEvidenceFactoryTest {

    private static final int MIN_BASE_QUALITY = 13;
    private final VariantHotspotEvidenceFactory victim = new VariantHotspotEvidenceFactory(0, MIN_BASE_QUALITY, Collections.emptySet());

    private static final String MNV_REF_SEQUENCE = "GATACAA";
    private static final VariantHotspot MNV =
            ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref("TAC").alt("CAT").build();

    private static final String SNV_REF_SEQUENCE = "GATAC";
    private static final VariantHotspot SNV =
            ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref("T").alt("C").build();
    private static final VariantHotspot INS =
            ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref("T").alt("TTT").build();
    private static final VariantHotspot DEL =
            ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref("TAC").alt("T").build();

    @Test
    public void testDelAlt() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfDelete(create(DEL), DEL, buildSamRecord(98, "3M2D1M", "GATA"));
        assertEvidence(evidence, 1, 1, 0, MIN_BASE_QUALITY, 1);
    }

    @Test
    public void testDelAlt2() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfDelete(create(DEL), DEL, buildSamRecord(98, "3M1D1M", "GATCA"));
        assertEvidence(evidence, 1, 0, 0, 0, 1);
    }

    @Test
    public void testDelAltAverageQuality() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfDelete(create(DEL),
                DEL,
                SAMRecordsTest.buildSamRecord(98, "3M2D1M", "GATA", buildQualities(13, 13, 16, 18)));
        assertEvidence(evidence, 1, 1, 0, 17, 1);
    }

    @Test
    public void testDelAltInsufficientQuality() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfDelete(create(DEL),
                DEL,
                SAMRecordsTest.buildSamRecord(98, "3M2D1M", "GATA", buildQualities(13, 13, 12, 12)));
        assertEvidence(evidence, 0, 0, 0, 0, 0);
    }

    @Test
    public void testDelAltNoBarrier() {
        final VariantHotspotEvidence evidence =
                victim.findEvidenceOfDelete(create(DEL), DEL, SAMRecordsTest.buildSamRecord(98, "3M2D", "GAT", buildQualities(13, 13, 16)));
        assertEvidence(evidence, 1, 1, 0, 16, 1);
    }

    @Test
    public void testDelIsActuallyInsert() {
        final VariantHotspotEvidence evidence =
                victim.findEvidenceOfDelete(create(DEL), DEL, SAMRecordsTest.buildSamRecord(98, "3M2I2M", "GATTTAC"));
        assertEvidence(evidence, 1, 0, 0, 0, 1);
    }

    @Test
    public void testDelIsRef() {
        final VariantHotspotEvidence evidence =
                victim.findEvidenceOfDelete(create(DEL), DEL, SAMRecordsTest.buildSamRecord(98, "5M", "GATAC"));
        assertEvidence(evidence, 1, 0, 1, 0, 0);
    }

    @Test
    public void testInsAlt() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "3M2I2M", "GATTTAC"));
        assertEvidence(evidence, 1, 1, 0, MIN_BASE_QUALITY, 1);
    }

    @Test
    public void testInsAlt2() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "3M2I2M", "GATATAC"));
        assertEvidence(evidence, 1, 0, 0, 0, 1);
    }

    @Test
    public void testInsAltAverageQuality() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS),
                INS,
                SAMRecordsTest.buildSamRecord(98, "3M2I2M", "GATTTAC", buildQualities(13, 13, 20, 9, 16, 13, 13)));
        assertEvidence(evidence, 1, 1, 0, 15, 1);
    }

    @Test
    public void testInsAltInsufficientQuality() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS),
                INS,
                SAMRecordsTest.buildSamRecord(98, "3M2I2M", "GATTTAC", buildQualities(13, 13, 12, 12, 12, 13, 13)));
        assertEvidence(evidence, 0, 0, 0, 0, 0);
    }

    @Test
    public void testInsAltNoBarrier() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "3M2I", "GATTT"));
        assertEvidence(evidence, 1, 1, 0, MIN_BASE_QUALITY, 1);
    }

    @Test
    public void testInsIsActuallyDel() {
        VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "3M2D1M", "GATA"));
        assertEvidence(evidence, 1, 0, 0, 0);

        evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "3M2D", "GAT"));
        assertEvidence(evidence, 1, 0, 0, 0, 1);
    }

    @Test
    public void testInsIsRef() {
        final VariantHotspotEvidence evidence = victim.findEvidenceOfInsert(create(INS), INS, buildSamRecord(98, "5M", "GATAC"));
        assertEvidence(evidence, 1, 0, 1, 0, 0);
    }

    @Test
    public void testSnvAlt() {
        final VariantHotspotEvidence evidence = createSNVEvidence(buildSamRecord(98, "5M", "GACAC"));
        assertEvidence(evidence, 1, 1, 0, 13);
    }

    @Test
    public void testMnvAlt() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAA"));
        assertEvidence(evidence, 1, 1, 0, 13);
    }

    @Test
    public void testMnvInsufficientQuality() {
        final VariantHotspotEvidence evidence =
                createMNVEvidence(SAMRecordsTest.buildSamRecord(98, "7M", "GACATAA", buildQualities(12, 12, 12, 12, 12, 12, 12)));
        assertEvidence(evidence, 0, 0, 0, 0);
    }

    @Test
    public void testMnvBarelySufficientQuality() {
        final VariantHotspotEvidence evidence =
                createMNVEvidence(SAMRecordsTest.buildSamRecord(98, "7M", "GACATAA", buildQualities(13, 13, 9, 20, 16, 13, 13)));
        assertEvidence(evidence, 1, 1, 0, 15);
    }

    @Test
    public void testMnvRef() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", MNV_REF_SEQUENCE));
        assertEvidence(evidence, 1, 0, 1, 0);
    }

    @Test
    public void testDelInMnv() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "3M1D3M", "GACTAA"));
        assertEvidence(evidence, 1, 0, 0, 0);
    }

    @Test
    public void testInsInMnv() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "4M1I3M", "GACATTAA"));
        assertEvidence(evidence, 1, 0, 0, 0);
    }

    @Test
    public void testInsInMnvButRefHasInsufficientQuality() {
        final VariantHotspotEvidence evidence =
                createMNVEvidence(SAMRecordsTest.buildSamRecord(98, "4M1I3M", "GACATTAA", buildQualities(13, 13, 12, 13, 13, 13, 13, 13)));
        assertEvidence(evidence, 0, 0, 0, 0);
    }

    @Test
    public void testInsBeforeMnv() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "2M1I5M", "GACCATAA"));
        assertEvidence(evidence, 1, 1, 0, 13);
    }

    @Test
    public void testMnvPartOfLarger() {
        VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GGCATAA"));
        assertEvidence(evidence, 1, 0, 0, 0);

        evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAC"));
        assertEvidence(evidence, 1, 0, 0, 0);
    }

    private VariantHotspotEvidence createMNVEvidence(SAMRecord record) {
        return victim.findEvidenceOfMNV(create(MNV), 98, MNV_REF_SEQUENCE, MNV, record);
    }

    private VariantHotspotEvidence createSNVEvidence(SAMRecord record) {
        return victim.findEvidenceOfMNV(create(SNV), 98, SNV_REF_SEQUENCE, SNV, record);
    }

    @Test
    public void testIsPartOfLargerMNV() {
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "GATACAA")));
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "GACATAA")));
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(97, "7M", "AGATACAAA")));
        assertFalse(isVariantPartOfLargerMNV(99, "ATACA", MNV, buildSamRecord(98, "7M", "TATACAG")));
        assertFalse(isVariantPartOfLargerMNV(99, "ATACA", MNV, buildSamRecord(98, "7M", "TATACAG")));

        // Indel or del in front
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "1M1D5M", "TTACAA")));
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "2M1I5M", "GGGTACAA")));

        // Indel or del behind
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "5M1I2M", "GATACCCC")));
        assertFalse(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "5M1D1M", "GATACC")));

        // Anything in the buffer is different
        assertTrue(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "AATACAA")));
        assertTrue(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "GGTACAA")));
        assertTrue(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "GATACAG")));
        assertTrue(isVariantPartOfLargerMNV(98, MNV_REF_SEQUENCE, MNV, buildSamRecord(98, "7M", "GATACGA")));
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString) {
        final StringBuilder qualityString = new StringBuilder();
        for (int i = 0; i < readString.length(); i++) {
            qualityString.append((char) (MIN_BASE_QUALITY + PHRED_OFFSET));
        }

        return SAMRecordsTest.buildSamRecord(alignmentStart, cigar, readString, qualityString.toString());
    }

    private static void assertEvidence(@NotNull VariantHotspotEvidence victim, int readDepth, int altCount, int refCount, int quality) {
        assertEquals(readDepth, victim.readDepth());
        assertEquals(altCount, victim.altSupport());
        assertEquals(refCount, victim.refSupport());
        assertEquals(quality, victim.altQuality());
    }

    private static void assertEvidence(@NotNull VariantHotspotEvidence victim, int readDepth, int altCount, int refCount, int quality,
            int indelSupport) {
        assertEvidence(victim, readDepth, altCount, refCount, quality);
        assertEquals(indelSupport, victim.indelSupport());
    }

    @NotNull
    private static String buildQualities(int... qualities) {
        final StringBuilder qualityString = new StringBuilder();
        for (final int quality : qualities) {
            qualityString.append((char) (quality + PHRED_OFFSET));
        }
        return qualityString.toString();
    }
}
