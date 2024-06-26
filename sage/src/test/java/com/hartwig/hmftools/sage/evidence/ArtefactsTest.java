package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.ArtefactContext;

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

        VariantReadContext readContext = createReadContext(
                variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        ArtefactContext artefactContext = ArtefactContext.buildContext(readContext);
        assertNotNull(artefactContext);

        assertFalse(artefactContext.hasHomopolymerOffset(SE_START));
        assertEquals(0, artefactContext.homopolymerOffset(SE_END));
        assertEquals((byte)'T', artefactContext.homopolymerBase(SE_END));

        int readLength = readContextBases.length();
        String cigar = format("%dM", readLength);
        byte[] readQualities = buildDefaultBaseQuals(readLength);

        byte lowQualBase = 10;
        readQualities[varIndex] = lowQualBase;
        SAMRecord record = buildSamRecord(90, cigar, readContextBases, new String(readQualities));

        byte adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(adjustedBaseQual, INVALID_BASE_QUAL);

        record.setReadNegativeStrandFlag(true);

        adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(lowQualBase, adjustedBaseQual);

        // single base DEL on left of HP
        variant = new SimpleVariant(CHR_1, pos, "CA", "C");

        varIndex = 16;
        readContextBases = REF_BASES.substring(0, varIndex + 1) + REF_BASES.substring(varIndex + 2);

        readContext = createReadContext(
                variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        artefactContext = ArtefactContext.buildContext(readContext);

        assertNotNull(artefactContext);
        assertEquals(1, artefactContext.homopolymerOffset(SE_END));
        assertEquals((byte)'T', artefactContext.homopolymerBase(SE_END));

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

        readContext = createReadContext(variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        artefactContext = ArtefactContext.buildContext(readContext);

        assertNotNull(artefactContext);
        assertEquals(2, artefactContext.homopolymerOffset(SE_END));
        assertEquals((byte)'T', artefactContext.homopolymerBase(SE_END));

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

        readContext = createReadContext(variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        artefactContext = ArtefactContext.buildContext(readContext);

        assertNotNull(artefactContext);
        assertEquals(3, artefactContext.homopolymerOffset(SE_END));
        assertEquals((byte)'T', artefactContext.homopolymerBase(SE_END));

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

        readContext = createReadContext(variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        artefactContext = ArtefactContext.buildContext(readContext);

        assertNotNull(artefactContext);
        assertEquals(-1, artefactContext.homopolymerOffset(SE_END));
        assertEquals((byte)'T', artefactContext.homopolymerBase(SE_END));

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

        VariantReadContext readContext = createReadContext(
                variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        ArtefactContext artefactContext = ArtefactContext.buildContext(readContext);

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
        assertEquals(INVALID_BASE_QUAL, adjustedBaseQual);

        // single base DEL on left of HP
        variant = new SimpleVariant(CHR_1, pos, "CA", "C");

        varIndex = 16;
        readContextBases = REF_BASES.substring(0, varIndex + 1) + REF_BASES.substring(varIndex + 2);

        readContext = createReadContext(
                variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        artefactContext = ArtefactContext.buildContext(readContext);


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
    public void testSnvWithAdjacentInsert()
    {
        int pos = 115;
        String leftFlank = "AAAATCAGGA";
        String rightFlank = "AAAGGCAGAA";

        SimpleVariant variant = new SimpleVariant(CHR_1, pos, "G", "T");

        int varIndex = 15;
        String readContextBases = leftFlank + "CTGTTTTTTTTTTTTTTTTTTA" + rightFlank;

        VariantReadContext readContext = createReadContext(
                variant, readContextBases, varIndex - 2, varIndex, varIndex + 10);

        ArtefactContext artefactContext = ArtefactContext.buildContext(readContext);

        int readLength = readContextBases.length();
        String cigar = format("%dM", readLength);
        byte[] readQualities = buildDefaultBaseQuals(readLength);

        byte lowQualBase = 10;
        readQualities[varIndex-1] = lowQualBase;
        SAMRecord record = buildSamRecord(20, cigar, readContextBases, new String(readQualities));
        record.setReadNegativeStrandFlag(true);

        byte adjustedBaseQual = artefactContext.findApplicableBaseQual(record, varIndex);

        assertEquals(lowQualBase, adjustedBaseQual);
    }
}
