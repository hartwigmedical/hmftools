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
        String refBases = BUFFER_REF_BASES + "ATTTTCGTCGT";
        setRefBases(refBases);

        // delete of 1 base
        SimpleVariant variant = new SimpleVariant(CHR_1, 10, "AT", "A");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        String readBases = BUFFER_REF_BASES + "ATTTCGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        short[] tpValues = new short[readBases.length()];
        tpValues[10] = 1;
        tpValues[12 ] = 1;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 9);
        assertEquals(34, calcQual);

        tpValues = new short[readBases.length()];
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
        tpValues = new short[readBases.length()];
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
        tpValues = new short[readBases.length()];
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
    }

    @Test
    public void testHomopolymerDeletion()
    {
        //                                    0123456789
        String refBases = BUFFER_REF_BASES + "ATTTTCGCAGT" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // the whole HP must be deleted
        SimpleVariant variant = new SimpleVariant(CHR_1, 10, "ATTTT", "A");

        UltimaQualModel model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        String readBases = BUFFER_REF_BASES + "ACGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        short[] tpValues = new short[readBases.length()];
        t0Values[9] = 25;
        t0Values[10] = 30;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 9);
        assertEquals(30, calcQual);

        read.setReadNegativeStrandFlag(true);

        calcQual = model.calculateQual(read, 9);
        assertEquals(30, calcQual);

        // test out of cycle deletions
        variant = new SimpleVariant(CHR_1, 16, "GC", "G");

        model = mModelBuilder.buildContext(variant);
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTCGAGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new short[readBases.length()];
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
        short[] tpValues = new short[readBases.length()];

        tpValues[11] = 2;
        tpValues[12] = 2;
        baseQualities[11] = 30;
        baseQualities[2] = 30;

        tpValues[13] = 2;
        tpValues[15] = 2;
        baseQualities[13] = 25;
        baseQualities[15] = 25;

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 12);
        assertEquals(17, calcQual);

        // as before but with longer, lopsided transitional delete
        variant = new SimpleVariant(CHR_1, 11, "TTTTA", "T");

        model = mModelBuilder.buildContext(variant);
        assertEquals(HOMOPOLYMER_TRANSITION, model.type());

        //                              0123456
        readBases = BUFFER_REF_BASES + "ATAAAAG";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new short[readBases.length()];

        tpValues[11] = 3;
        baseQualities[11] = 30;

        tpValues[12] = -1;
        tpValues[13] = 1;
        tpValues[14] = 1;
        tpValues[15] = -1;
        baseQualities[13] = 35;
        baseQualities[14] = 35;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 11);
        assertEquals(27, calcQual);
    }

    private static SAMRecord buildUltimaRead(
            final String readBases, final int readStart, final byte[] qualities, final short[] tpValues, final byte[] t0Values)
    {
        String cigar = format("%dM", qualities.length);
        SAMRecord record = buildSamRecord(readStart, cigar, readBases, qualities);
        record.setAttribute(TP_TAG, tpValues);
        record.setAttribute(T0_TAG, new String(t0Values));

        return record;
    }
}
