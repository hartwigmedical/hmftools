package com.hartwig.hmftools.sage.quality;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.TP_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.HOMOPOLYMER_ADJUSTMENT;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.HOMOPOLYMER_DELETION;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.HOMOPOLYMER_TRANSITION;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.SNV;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelTest
{
    private final MockRefGenome mRefGenome;
    private final UltimaQualCalculator mModelBuilder;

    private static final String BUFFER_REF_BASES = "AACCGGTTA";

    public UltimaQualModelTest()
    {
        mRefGenome = new MockRefGenome(true);
        mModelBuilder = new UltimaQualCalculator(mRefGenome);
    }

    private void setRefBases(final String refBases)
    {
        mRefGenome.RefGenomeMap.put(CHR_1, refBases);
    }

    @Test
    public void testHomopolymerAdjustment()
    {
        //                                    0123456789
        String refBases = BUFFER_REF_BASES + "ATTTTCGTCGT" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // delete of 1 base
        SimpleVariant variant = new SimpleVariant(CHR_1, 10, "AT", "A");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        String readBases = BUFFER_REF_BASES + "ATTTCGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];
        tpValues[10] = 1;
        tpValues[12 ] = 1;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 9);
        assertEquals(34, calcQual);

        tpValues = new byte[readBases.length()];
        tpValues[11] = 1;
        baseQualities[11] = 25;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 9);
        assertEquals(25, calcQual);

        // delete of 3 bases
        variant = new SimpleVariant(CHR_1, 10, "ATTT", "A");

        model = mModelBuilder.buildContext(variant);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATCGTCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        tpValues[10] = 3;
        baseQualities[10] = 25;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 9);
        assertEquals(25, calcQual);

        // test where the required adjustment value isn't in the TP values
        tpValues[10] = -1;

        calcQual = model.calculateQual(read, 9);
        assertEquals(40, calcQual);

        // insert of 2 bases
        variant = new SimpleVariant(CHR_1, 10, "A", "ATT");

        model = mModelBuilder.buildContext(variant);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTTTCGTCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        tpValues[10] = 1;
        tpValues[11] = -1;
        tpValues[12] = -2; // will use these 2
        tpValues[13] = -2;
        baseQualities[12] = 30;
        baseQualities[13] = 30;
        tpValues[14] = 11;
        tpValues[15] = 1;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 9);
        assertEquals(27, calcQual);

        // C>CA in TT C GTC, insert of base which doesn't match the existing context
        variant = new SimpleVariant(CHR_1, 15, "C", "CA");

        model = mModelBuilder.buildContext(variant);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTCAGTCGT" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        tpValues[16] = -1;
        baseQualities[16] = 32;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 15);
        assertEquals(32, calcQual);
    }

    @Test
    public void testHomopolymerDeletion()
    {
        //                                    0123456789
        String refBases = BUFFER_REF_BASES + "ATTTTCGACGT" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // the whole HP must be deleted
        SimpleVariant variant = new SimpleVariant(CHR_1, 10, "ATTTT", "A");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        String readBases = BUFFER_REF_BASES + "ACGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];
        t0Values[9] = 25;
        t0Values[10] = 30;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 9);
        assertEquals(30, calcQual);

        read.setReadNegativeStrandFlag(true);

        calcQual = model.calculateQual(read, 9);
        assertEquals(30, calcQual);

        // test out of cycle deletions
        variant = new SimpleVariant(CHR_1, 16, "GA", "G");

        model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTCGCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 9);
        assertEquals(ULTIMA_MAX_QUAL, calcQual);

        read.setReadNegativeStrandFlag(true);
        calcQual = model.calculateQual(read, 9);
        assertEquals(ULTIMA_MAX_QUAL, calcQual);
    }

    @Test
    public void testHomopolymerTransition()
    {
        //                                    01234567890
        String refBases = BUFFER_REF_BASES + "ATTTTAAAAAG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // a deletion which crosses 2 HPs
        SimpleVariant variant = new SimpleVariant(CHR_1, 12, "TTTAA", "T");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_TRANSITION, model.type());

        //                                     0123456
        String readBases = BUFFER_REF_BASES + "ATTAAAG";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];

        tpValues[11] = 2;
        tpValues[12] = 2;
        baseQualities[11] = 11;
        baseQualities[2] = 11;

        tpValues[13] = 2;
        tpValues[15] = 2;
        baseQualities[13] = 25;
        baseQualities[15] = 25;

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 12);
        assertEquals(30, calcQual);

        // as before but with longer, lopsided transitional delete
        variant = new SimpleVariant(CHR_1, 11, "TTTTA", "T");

        model = mModelBuilder.buildContext(variant);
        assertEquals(HOMOPOLYMER_TRANSITION, model.type());

        //                              0123456
        readBases = BUFFER_REF_BASES + "ATAAAAG";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        tpValues[11] = 3;
        baseQualities[11] = 25;

        tpValues[12] = -1;
        tpValues[13] = 1;
        tpValues[14] = 1;
        tpValues[15] = -1;
        baseQualities[13] = 12;
        baseQualities[14] = 12;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 11);
        assertEquals(34, calcQual);
    }

    @Test
    public void testSNVs()
    {
        //                                    01     234     56
        String refBases = BUFFER_REF_BASES + "AG" + "ACA" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // C>T in ACA > ATA, matches neither side's ref so reverts to max qual
        SimpleVariant variant = new SimpleVariant(CHR_1, 13, "C", "T");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(SNV, model.type());

        String readBases = BUFFER_REF_BASES + "AGATAAG";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 13);
        assertEquals(37, calcQual);

        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "CCG" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // C>T in CCG > CTG, left contraction, right insertion/expansion
        variant = new SimpleVariant(CHR_1, 13, "C", "T");

        model = mModelBuilder.buildContext(variant);

        readBases = BUFFER_REF_BASES + "AGCTGAG";

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the deleted base qual
        tpValues[12] = 1;
        baseQualities[12] = 16;

        tpValues[13] = -1;
        baseQualities[13] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 13);
        assertEquals(37, calcQual);

        // C>T GCT > GTT, left full delete, right insertion/expansion
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "GCT" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant);

        //                              0123456
        readBases = BUFFER_REF_BASES + "AGGTTAG" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        t0Values[12] = 10;
        t0Values[13] = 15; // an SNV with any HP deletion, including of a single base, just uses t0

        tpValues[13] = -1;
        baseQualities[13] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 13);
        assertEquals(15, calcQual);

        // C>T in TCG > TTG, left 1-base ins/expansion, right full delete
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "TCG" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant);

        //                              0123456
        readBases = BUFFER_REF_BASES + "AGTTGAG" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        t0Values[13] = 25;
        t0Values[14] = 10;

        tpValues[12] = tpValues[13] = -1;
        baseQualities[12] = baseQualities[13] = 20;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 13);
        assertEquals(25, calcQual);


        // C>T in GCC > GTC, left contraction, right insertion/expansion
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "GCC" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant);

        readBases = BUFFER_REF_BASES + "AGGTCAG";

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the delete base qual
        tpValues[14] = 1;
        baseQualities[14] = 16;

        tpValues[13] = -1;
        baseQualities[13] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 13);
        assertEquals(37, calcQual);


        // contraction on the C side, expansion on the T side
        // C>T in TTT TCC CCC > TTT TTC CCC

        variant = new SimpleVariant(CHR_1, 14, "C", "T");

        //                             012     345     678
        refBases = BUFFER_REF_BASES + "TTT" + "TCC" + "CCC" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant);

        //                              012345678
        readBases = BUFFER_REF_BASES + "TTTTTCCCC" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the HP insert base quals
        tpValues[10] = tpValues[14] = -1;
        baseQualities[10] = baseQualities[14] = 16;

        // the HP delete quals
        tpValues[15] = tpValues[18] = 1;
        baseQualities[15] = baseQualities[18] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 14);
        assertEquals(31, calcQual);

        // HP on right partially deleted -> right matches ref, insert on left
        // AG GCC AGA > AG GTC AGA

        // HP on at var position fully deleted -> right matches alt, like a HP expansion on right
        // AG GCT AGA > AG GTT AGA

        // HP at var position fully deleted -> left matches alt, like a HP expansion on left
        // AG TCG AGA > AG TTG AGA


        // variant has one straddling base matching ref or alt
        // AG CCG AGA > AG CTG AGA
        // delete of C on left,  insert of T


        // variant has one straddling base matching ref or alt
        // TC TCG GC > TC TTG GC
        // full delete of C, insert of additional T on left
    }

    private static SAMRecord buildUltimaRead(
            final String readBases, final int readStart, final byte[] qualities, final byte[] tpValues, final byte[] t0Values)
    {
        String cigar = format("%dM", qualities.length);
        SAMRecord record = buildSamRecord(readStart, cigar, readBases, qualities);
        record.setAttribute(TP_TAG, tpValues);
        record.setAttribute(T0_TAG, new String(t0Values));

        return record;
    }
}
