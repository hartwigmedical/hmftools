package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.consensus.UltimaRoutines.formLowQualTag;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.test.SeqTechTestUtils;

import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaReadTest
{
    private static final byte LOW_QUAL = 11;

    @Ignore
    @Test
    public void testUltimaQualTag()
    {
        int readStart = 100;
        String readBases = "ATTTC";

        byte[] baseQuals = buildBaseQuals(readBases.length(), ULTIMA_MAX_QUAL);
        byte[] tpValues = new byte[readBases.length()];
        byte[] t0Values = buildBaseQuals(readBases.length(), ULTIMA_MAX_QUAL);

        // test 0: no quals low enough to trigger creation of the tag

        SAMRecord read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        String ulqTag = formLowQualTag(read);

        assertNull(ulqTag);

        byte lowQualBase = 12;

        // test 1:
        // bases: A  T  T  T  C
        // index: 0  1  2  3  4
        // qual:  35 12 35 12 35
        // t0:    35 35 35 35 35
        // tp:    0  -1 0 -1  0

        resetQualValues(baseQuals, tpValues, t0Values);

        setQualValues(baseQuals, List.of(1, 3), lowQualBase);

        read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        ulqTag = formLowQualTag(read);

        assertNotNull(ulqTag);
        assertEquals("1-3", ulqTag);

        // test 2:
        // bases: A  T  T  T  C
        // index: 0  1  2  3  4
        // qual:  35 12 35 12 35
        // t0:    35 10 10 10 10
        // tp:    0  -1 0 -1  0

        resetQualValues(baseQuals, tpValues, t0Values);

        setQualValues(baseQuals, List.of(1, 3), lowQualBase);
        setQualValues(t0Values, List.of(1, 2, 3, 4), lowQualBase);

        read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        ulqTag = formLowQualTag(read);

        assertNotNull(ulqTag);
        assertEquals("1-4", ulqTag);

        // test 3:
        // bases: A  G  T  T  T  C
        // index: 0  1  2  3  4  5
        // qual:  35 35 23 35 23 23
        // t0:    12 12 35 35 35 35
        // tp:    0  0  -1 0 -1  0

        readBases = "AGTTTC";

        baseQuals = buildBaseQuals(readBases.length(), ULTIMA_MAX_QUAL);
        tpValues = new byte[readBases.length()];
        t0Values = buildBaseQuals(readBases.length(), ULTIMA_MAX_QUAL);

        setQualValues(baseQuals, List.of(2, 4, 5), (byte)23);
        setQualValues(t0Values, List.of(0, 1), lowQualBase);

        read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        ulqTag = formLowQualTag(read);

        assertNotNull(ulqTag);
        assertEquals("0-1,5", ulqTag);
    }

    private static void resetQualValues(final byte[] quals, final byte[] tpValues, final byte[] t0Values)
    {
        for(int i = 0; i < quals.length; ++i)
        {
            quals[i] = ULTIMA_MAX_QUAL;
            tpValues[i] = 0;
            t0Values[i] = ULTIMA_MAX_QUAL;
        }
    }

    private static void setQualValues(final byte[] quals, final List<Integer> indices, final byte qualValue)
    {
        for(Integer index : indices)
        {
            quals[index] = qualValue;
        }
    }

}
