package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.sage.common.TestUtils.QUALITY_CALCULATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.RECALIBRATION;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.evidence.ArtefactContext.NOT_APPLICABLE_BASE_QUAL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ArtefactsTest
{
    private static final String FLANK_BASES = "AAAAAGGGGG";
    private static final String HOMOPOLYMER_SEQ = "TTTTTTTT";

    private static final String REF_BASES = FLANK_BASES + "ACGTTGCA" + HOMOPOLYMER_SEQ + "ACGTTGCA" + FLANK_BASES;
    //            10         20         30
    // 0123456789 01234567 89012345 67890123
    // AAAAAGGGGG ACGTTGCA TTTTTTTT ACGTTGCA AAAAAGGGGG

    @Test
    public void testHomopolymerArtefacts()
    {
        int pos = 100;

        // A>T where SNV is immediately next to (left of) upstream homopolymer (HP)
        SimpleVariant variant = new SimpleVariant(CHR_1, pos, "A", "T");

        int varIndex = 17;
        String readContextBases = REF_BASES.substring(0, varIndex) + "T" + REF_BASES.substring(varIndex + 1);

        IndexedBases indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        ReadContext readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        ReadContextCounter rcCounter = new ReadContextCounter(
                1, variant, readContext, VariantTier.PANEL, 100, 0,
                TEST_CONFIG, QUALITY_CALCULATOR, null);

        ArtefactContext artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);

        assertFalse(artefactContext.hasHomopolymerOffset(SE_START));
        assertEquals(0, artefactContext.homopolymerOffset(SE_END));

        int readLength = readContextBases.length();
        String cigar = format("%dM", readLength);
        byte[] readQualities = buildDefaultBaseQuals(readLength);

        byte lowQualBase = 10;
        readQualities[varIndex] = lowQualBase;
        SAMRecord record = buildSamRecord(90, cigar, readContextBases, new String(readQualities));

        byte adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(adjustedBaseQual, NOT_APPLICABLE_BASE_QUAL);

        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(lowQualBase, adjustedBaseQual);

        // single base DEL on left of HP
        variant = new SimpleVariant(CHR_1, pos, "CA", "C");

        varIndex = 16;
        readContextBases = REF_BASES.substring(0, varIndex + 1) + REF_BASES.substring(varIndex + 2);

        indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertEquals(1, artefactContext.homopolymerOffset(SE_END));

        readQualities = buildDefaultBaseQuals(readLength);
        readQualities[varIndex + 1] = lowQualBase;
        record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(lowQualBase, adjustedBaseQual);

        // similar before but where the variant is left aligned so not immediately next to the HP
        variant = new SimpleVariant(CHR_1, pos, "TG", "T");

        varIndex = 15;
        //             0123456 7890
        // ACGTTTGG -> ACGTTTG TTTT

        //            10             20         30
        // 0123456789 012345     678901234 567890123
        // AAAAAGGGGG ACGTTT [G] GTTTTTTTT ACGTTGCA AAAAAGGGGG
        // low-qual               X
        readContextBases = FLANK_BASES + "ACGTTTG" + REF_BASES.substring(18);

        indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertEquals(2, artefactContext.homopolymerOffset(SE_END));

        readQualities = buildDefaultBaseQuals(readLength);
        readQualities[17] = lowQualBase;
        record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(lowQualBase, adjustedBaseQual);

        // similar to above but the delete is 2-bases left-aligned, so not immediately left of the HP
        //            10              20         30
        // 0123456789 01234567      89012345 67890123
        // AAAAAGGGGG AACGTAGT [GT] TTTTTTTT ACGTTGCA AAAAAGGGGG
        // low-qual                 X
        variant = new SimpleVariant(CHR_1, pos, "AGT", "A");

        varIndex = 15;
        readContextBases = FLANK_BASES + "AACGTAGT" + HOMOPOLYMER_SEQ + REF_BASES.substring(26);

        indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertEquals(3, artefactContext.homopolymerOffset(SE_END));

        readQualities = buildDefaultBaseQuals(readLength);
        readQualities[18] = lowQualBase;
        record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(lowQualBase, adjustedBaseQual);


        // with germline insertion which needs to be ignored when looking for the ref base
        // AAAAAGGGGG ACGTTGC insT [T] TTTTTTTT ACGTTGCA AAAAAGGGGG   - read
        // AAAAAGGGGG ACGTTGC      [A] TTTTTTTT ACGTTGCA AAAAAGGGGG   - ref

        // an SNV which follows a germline T insertion, which should be ignored when finding the first ref base
        //            10                20         30
        // 0123456789 0123456    7  8  9  01234567 890123
        // AAAAAGGGGG AACGTAG insT [T] TTTTTTTT ACGTTGCA AAAAAGGGGG
        // low-qual                 X
        variant = new SimpleVariant(CHR_1, pos, "A", "T");

        varIndex = 18;
        readContextBases = FLANK_BASES + "AACGTAG" + "T" + "T" + HOMOPOLYMER_SEQ + REF_BASES.substring(26);

        indexBases = new IndexedBases(
                pos, varIndex, 15, 27, 10, readContextBases.getBytes());

        readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertEquals(-1, artefactContext.homopolymerOffset(SE_END));

        readQualities = buildDefaultBaseQuals(readLength);
        readQualities[18] = lowQualBase;
        cigar = "17M1I27M";
        record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(lowQualBase, adjustedBaseQual);
    }

    @Test
    public void testHomopolymerArtefactsPosStrand()
    {
        int pos = 100;

        SimpleVariant variant = new SimpleVariant(CHR_1, pos, "A", "T");

        int varIndex = 26;
        String readContextBases = REF_BASES.substring(0, varIndex) + "T" + REF_BASES.substring(varIndex + 1);

        IndexedBases indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        ReadContext readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        ReadContextCounter rcCounter = new ReadContextCounter(
                1, variant, readContext, VariantTier.PANEL, 100, 0,
                TEST_CONFIG, QUALITY_CALCULATOR, null);

        ArtefactContext artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertFalse(artefactContext.hasHomopolymerOffset(SE_END));
        assertEquals(0, artefactContext.homopolymerOffset(SE_START));

        int readLength = readContextBases.length();
        String cigar = format("%dM", readLength);
        byte[] readQualities = buildDefaultBaseQuals(readLength);

        byte lowQualBase = 10;
        readQualities[varIndex] = lowQualBase;
        SAMRecord record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));

        byte adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(lowQualBase, adjustedBaseQual);

        // ignore -ve strand for downstream HPs
        record.setReadNegativeStrandFlag(true);
        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(NOT_APPLICABLE_BASE_QUAL, adjustedBaseQual);

        // single base DEL on left of HP
        variant = new SimpleVariant(CHR_1, pos, "CA", "C");

        varIndex = 16;
        readContextBases = REF_BASES.substring(0, varIndex + 1) + REF_BASES.substring(varIndex + 2);

        indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 27, 10, readContextBases.getBytes());

        readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        artefactContext = ArtefactContext.buildContext(variant, indexBases);
        assertNotNull(artefactContext);
        assertEquals(1, artefactContext.homopolymerOffset(SE_END));

        readQualities = buildDefaultBaseQuals(readLength);
        readQualities[varIndex + 1] = lowQualBase;
        record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);
        assertEquals(lowQualBase, adjustedBaseQual);
    }

    @Test
    public void testSnvBeforeInsert()
    {
        // create an SNV directly before a 1-base insert - confirm impact on read contexts and how reads are handled
        int pos = 117;

        String refBases = FLANK_BASES + "ACGTTCCAACCTTGCA" + FLANK_BASES;
        //            10                20          30
        // 0123456789 0123456 7      8 9012345  6789012345
        // AAAAAGGGGG ACGTTCC A>G insT ACCTTGCA AAAAAGGGGG

        SimpleVariant variant = new SimpleVariant(CHR_1, pos, "A", "G");

        int varIndex = 17;
        String readContextBases = refBases.substring(0, varIndex) + "G" + "T" + refBases.substring(varIndex + 1);

        IndexedBases indexBases = new IndexedBases(
                pos, varIndex, varIndex - 2, 22, FLANK_BASES.length(), readContextBases.getBytes());

        ReadContext readContext = new ReadContext(pos, "", 0, "", indexBases, false);

        IndexedBases indexedRefBases = new IndexedBases(100, 0, refBases.getBytes());
        QualityCalculator qualityCalculator = new QualityCalculator(TEST_CONFIG, RECALIBRATION, indexedRefBases);

        ReadContextCounter rcCounter = new ReadContextCounter(
                1, variant, readContext, VariantTier.PANEL, 100, 0,
                TEST_CONFIG, qualityCalculator, null);

        String readBases = readContextBases;
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String cigar = format("%dM",readBases.length());

        // the read's alignment start with the first base of the read context
        SAMRecord record = buildSamRecord(100, cigar, readContextBases, new String(baseQuals));

        rcCounter.processRead(record, 1, null);

        assertEquals(1, rcCounter.altSupport());

        cigar = "18M1I17M";

        record = buildSamRecord(100, cigar, readContextBases, new String(baseQuals));

        rcCounter.processRead(record, 1, null);

        // TODO: address this in CigarTraversal
        // assertEquals(1, rcCounter.altSupport());
    }
}
