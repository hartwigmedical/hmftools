package com.hartwig.hmftools.sage.quality;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.TP_TAG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.quality.UltimaVariantModel.HOMOPOLYMER_ADJUSTMENT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelTest
{
    private final MockRefGenome mRefGenome;
    private final UltimaQualModelBuilder mModelBuilder;

    private static final String BUFFER_REF_BASES = "AACCGGTTA";

    public UltimaQualModelTest()
    {
        mRefGenome = new MockRefGenome(true);
        mModelBuilder = new UltimaQualModelBuilder(mRefGenome);
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
        tpValues[13] = 1;
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

    private static SAMRecord buildUltimaRead(
            final String readBases, final int readStart, final byte[] qualities, final short[] tpValues, final byte[] t0Values)
    {
        String cigar = format("%dM", qualities.length);
        SAMRecord record = buildSamRecord(readStart, cigar, readBases, qualities);
        record.setAttribute(TP_TAG, tpValues);
        record.setAttribute(T0_TAG, t0Values);

        return record;
    }
}
