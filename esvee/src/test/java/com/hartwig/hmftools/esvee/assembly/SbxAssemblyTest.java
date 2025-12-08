package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CIGAR_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setSbxSequencing;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.passSbxDistinctPrimePositionsFilter;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.trimSbxUncertainBases;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHigherBaseQualCategory;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

import org.junit.After;
import org.junit.Test;

public class SbxAssemblyTest
{
    public SbxAssemblyTest()
    {
        setSbxSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testBaseQualComparisons()
    {
        byte qual1 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        byte qual2 = BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD;
        assertTrue(isHigherBaseQualCategory(qual1, qual2));

        qual1 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        qual2 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        assertFalse(isHigherBaseQualCategory(qual1, qual2));

        qual2 = SbxBamUtils.SBX_DUPLEX_QUAL;
        assertFalse(isHigherBaseQualCategory(qual1, qual2));
        assertTrue(isHigherBaseQualCategory(qual2, qual1));

        // assertFalse(isHigherBaseQualCategory(qual1, qual2));
    }

    @Test
    public void testTrimUncertainBases()
    {
        String readBases = REF_BASES_RANDOM_100;

        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());

        // low qual will be 70% of outer bases
        for(int i = 3; i < 10; ++i)
        {
            baseQualities[i] = BaseQualAdjustment.BASE_QUAL_MINIMUM;
            baseQualities[baseQualities.length - i - 1] = BaseQualAdjustment.BASE_QUAL_MINIMUM;
        }

        Read read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        trimSbxUncertainBases(read);

        assertEquals(10, read.trimCountStart());
        assertEquals(10, read.trimCountEnd());

        // now with a read that can trim from both sides past each other
        baseQualities = buildDefaultBaseQuals(readBases.length());

        for(int i = 0; i < 80; ++i)
        {
            baseQualities[i] = BaseQualAdjustment.BASE_QUAL_MINIMUM;
        }

        read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        trimSbxUncertainBases(read);

        assertEquals(80, read.trimCountStart());
        assertEquals(0, read.trimCountEnd());

        // and the other side
        baseQualities = buildDefaultBaseQuals(readBases.length());

        for(int i = 20; i < 100; ++i)
        {
            baseQualities[i] = BaseQualAdjustment.BASE_QUAL_MINIMUM;
        }

        read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        trimSbxUncertainBases(read);

        assertEquals(0, read.trimCountStart());
        assertEquals(80, read.trimCountEnd());
    }

    @Test
    public void testReadPrimeRange()
    {
        String readBases = REF_BASES_RANDOM_100;

        String cigar = "60M40S";

        List<Read> reads = Lists.newArrayList(
                createRead(TEST_READ_ID, 100, readBases, cigar),
                createRead(TEST_READ_ID, 101, readBases, cigar),
                createRead(TEST_READ_ID, 102, readBases, cigar));

        int juncReadDistance = 60; // irrelevant for this test
        int matchMismatch = 0;

        List<SupportRead> supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, juncReadDistance, matchMismatch, matchMismatch))
                .collect(Collectors.toList());

        assertFalse(passSbxDistinctPrimePositionsFilter(supportReads));

        Read read = createRead(TEST_READ_ID, 104, readBases, cigar);
        supportReads.add(new SupportRead(read, SupportType.JUNCTION, juncReadDistance, matchMismatch, matchMismatch ));

        assertTrue(passSbxDistinctPrimePositionsFilter(supportReads));

        // test reverse orientation reads
        cigar = "40S60M";
        reads = Lists.newArrayList(
                createRead(TEST_READ_ID, 102, readBases, cigar),
                createRead(TEST_READ_ID, 101, readBases, cigar),
                createRead(TEST_READ_ID, 100, readBases, cigar));

        reads.forEach(x -> x.bamRecord().setReadNegativeStrandFlag(true));

        supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, juncReadDistance, matchMismatch, matchMismatch))
                .collect(Collectors.toList());

        assertFalse(passSbxDistinctPrimePositionsFilter(supportReads));

        read = createRead(TEST_READ_ID, 98, readBases, cigar);
        read.bamRecord().setReadNegativeStrandFlag(true);
        supportReads.add(new SupportRead(read, SupportType.JUNCTION, juncReadDistance, matchMismatch, matchMismatch ));

        assertTrue(passSbxDistinctPrimePositionsFilter(supportReads));
    }
}
