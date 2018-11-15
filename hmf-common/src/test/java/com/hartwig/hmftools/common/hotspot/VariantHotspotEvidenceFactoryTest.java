package com.hartwig.hmftools.common.hotspot;

import static com.hartwig.hmftools.common.hotspot.VariantHotspotEvidenceFactory.isVariantPartOfLargerMNV;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantHotspotEvidenceFactoryTest {

    private static final String MNV_REF_SEQUENCE = "GATACAA";
    private static final VariantHotspot MNV = ImmutableVariantHotspot.builder().chromosome("11").position(100).ref("TAC").alt("CAT").build();

    private static final String SNV_REF_SEQUENCE = "GATAC";
    private static final VariantHotspot SNV = ImmutableVariantHotspot.builder().chromosome("11").position(100).ref("T").alt("C").build();

    @Test
    public void testSnvAlt() {
        final VariantHotspotEvidence evidence = createSNVEvidence(buildSamRecord(98, "7M", "GACAC"));
        assertEvidence(evidence, 1, 1, 0, 13);
    }

    @Test
    public void testMnvAlt() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAA"));
        assertEvidence(evidence, 1, 1, 0, 39);
    }

    @Test
    public void testMnvInsufficientQuality() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAA", buildQualities(12, 12, 12, 12, 12, 12, 12)));
        assertEvidence(evidence, 0, 0, 0, 0);
    }

    @Test
    public void testMnvBarelySufficientQuality() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAA", buildQualities(13, 13, 12, 20, 20, 13, 13)));
        assertEvidence(evidence, 1, 1, 0, 52);
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
    public void testInsInMnvButRefHasInsufficentQuality() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "4M1I3M", "GACATTAA", buildQualities(13, 13, 12, 13, 13, 13, 13, 13)));
        assertEvidence(evidence, 0, 0, 0, 0);
    }

    @Test
    public void testInsBeforeMnv() {
        final VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "2M1I5M", "GACCATAA"));
        assertEvidence(evidence, 1, 1, 0, 39);
    }

    @Test
    public void testMnvPartOfLarger() {
        VariantHotspotEvidence evidence = createMNVEvidence(buildSamRecord(98, "7M", "GGCATAA"));
        assertEvidence(evidence, 1, 0, 0, 0);

        evidence = createMNVEvidence(buildSamRecord(98, "7M", "GACATAC"));
        assertEvidence(evidence, 1, 0, 0, 0);
    }

    private VariantHotspotEvidence createMNVEvidence(SAMRecord record) {
        return VariantHotspotEvidenceFactory.findEvidenceOfMNV(VariantHotspotEvidenceFactory.create(MNV), 98, MNV_REF_SEQUENCE, MNV, record);
    }

    private VariantHotspotEvidence createSNVEvidence(SAMRecord record) {
        return VariantHotspotEvidenceFactory.findEvidenceOfMNV(VariantHotspotEvidenceFactory.create(SNV), 98, SNV_REF_SEQUENCE, SNV, record);
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

    @Test
    public void testInsertedBasesAfterPosition() {
        assertEquals(0, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(1, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "3M1I3M", "GATTACA")));
        assertEquals(2, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "3M2I3M", "GATTTACA")));

        assertEquals(0, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(0, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "3M1D2M", "GATCA")));
        assertEquals(0, VariantHotspotEvidenceFactory.insertedBasesAfterPosition(100, buildSamRecord(98, "3M2D1M", "GATA")));
    }

    @Test
    public void testDeletedBasesAfterPosition() {
        assertEquals(0, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(0, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "3M1I3M", "GATTACA")));
        assertEquals(0, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "3M2I3M", "GATTTACA")));

        assertEquals(0, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(1, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "3M1D2M", "GATCA")));
        assertEquals(2, VariantHotspotEvidenceFactory.deletedBasesAfterPosition(100, buildSamRecord(98, "3M2D1M", "GATA")));
    }

    @NotNull
    private static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString) {
        final StringBuilder qualityString = new StringBuilder();
        for (int i = 0; i < readString.length(); i++) {
            qualityString.append((char) (VariantHotspotEvidenceFactory.MIN_BASE_QUALITY + VariantHotspotEvidenceFactory.PHRED_OFFSET));
        }

        return buildSamRecord(alignmentStart, cigar, readString, qualityString.toString());
    }

    @NotNull
    private static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities) {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }

    private static void assertEvidence(@NotNull VariantHotspotEvidence victim, int readDepth, int altCount, int refCount, int quality) {
        assertEquals(readDepth, victim.readDepth());
        assertEquals(altCount, victim.altSupport());
        assertEquals(refCount, victim.refSupport());
        assertEquals(quality, victim.quality());
    }

    private static String buildQualities(int... qualities) {
        final StringBuilder qualityString = new StringBuilder();
        for (final int quality : qualities) {
            qualityString.append((char) (quality + VariantHotspotEvidenceFactory.PHRED_OFFSET));
        }
        return qualityString.toString();
    }

}
