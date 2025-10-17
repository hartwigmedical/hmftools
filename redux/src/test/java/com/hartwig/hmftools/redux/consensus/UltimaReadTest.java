package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG_DELIM;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.consensus.UltimaRoutines.setLowQualTag;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

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
        String readBases = "AACCGGTTACGT";

        byte[] baseQuals = buildBaseQuals(readBases.length(), ULTIMA_MAX_QUAL);
        byte[] tpValues = new byte[readBases.length()];
        byte[] t0Values = new byte[readBases.length()];

        // test 1: no quals low enough to trigger creation of the tag
        SAMRecord read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        setLowQualTag(read);

        String ulqTag = read.getStringAttribute(ULT_QUAL_TAG_DELIM);

        assertNull(ulqTag);

        // qual designations: L = low, M = medium, H = max qual


        // test 2: complete deletion only
        // index: 012345678901
        // t0:    HHLLHHHLLHHH
        // tp:    all max
        // qual:

        t0Values[2] = LOW_QUAL;
        t0Values[3] = LOW_QUAL;

        t0Values[7] = LOW_QUAL;
        t0Values[8] = LOW_QUAL;

        read = SeqTechTestUtils.buildUltimaRead(READ_ID_GEN.nextId(), CHR_1, readStart, readBases, baseQuals, tpValues, t0Values);

        setLowQualTag(read);

        ulqTag = read.getStringAttribute(ULT_QUAL_TAG_DELIM);

        assertNotNull(ulqTag);

    }

}
