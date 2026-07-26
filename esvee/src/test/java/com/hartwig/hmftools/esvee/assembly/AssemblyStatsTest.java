package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.JunctionReadTypes.calcProximateJunctionReadRatio;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.junit.Test;

public class AssemblyStatsTest
{
    @Test
    public void testAssemblySupplementaryStats()
    {
        Junction negJunction = new Junction(CHR_1, 100, REVERSE);

        String assemblyBases = REF_BASES_200.substring(0, 100);
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly assembly = new JunctionAssembly(negJunction, assemblyBases.getBytes(), baseQuals, 50);

        int readStartPos = 100;
        Read read = createRead(TEST_READ_ID, readStartPos, assemblyBases, "40S60M");
        read.bamRecord().setAttribute(SUPPLEMENTARY_ATTRIBUTE, "1,1000,+,40M60S,30,0");

        assembly.addJunctionRead(read);

        read = createRead(TEST_READ_ID, readStartPos, assemblyBases, "40S60M");
        read.bamRecord().setAttribute(SUPPLEMENTARY_ATTRIBUTE, "1,1400,+,40M60S,30,0");

        assembly.addJunctionRead(read);

        read = createRead(TEST_READ_ID, readStartPos, assemblyBases, "40S60M");
        read.bamRecord().setAttribute(SUPPLEMENTARY_ATTRIBUTE, "1,900,+,40M60S,25,0"); // below min alignment score

        assembly.addJunctionRead(read);

        read = createRead(TEST_READ_ID, readStartPos, assemblyBases, "40S60M");
        read.bamRecord().setAttribute(SUPPLEMENTARY_ATTRIBUTE, "2,1000,+,40M60S,30,0");

        assembly.addJunctionRead(read);

        for(SupportRead supportRead : assembly.support())
        {
            assembly.stats().addRead(supportRead, negJunction, supportRead.cachedRead());
        }

        assertEquals(0.67, assembly.stats().suppRemoteRegionRatio(), 0.1);
    }

    @Test
    public void testAssemblyProx()
    {
        Junction negJunction = new Junction(CHR_1, 100, REVERSE);

        String assemblyBases = REF_BASES_200.substring(0, 100);

        List<Read> rawReads = Lists.newArrayList(
                createRead(TEST_READ_ID, 100, assemblyBases, "40S60M"), // exact and long enough
                createRead(TEST_READ_ID, 100, assemblyBases, "35S65M"), // as above
                createRead(TEST_READ_ID, 100, assemblyBases, "20S80M"), // too short
                createRead(TEST_READ_ID, 70, assemblyBases, "10S90M"), // too far away to be counted
                createRead(TEST_READ_ID, 90, assemblyBases, "40S60M")); // included but not exact

        assertEquals(0.5, calcProximateJunctionReadRatio(negJunction, rawReads), 0.1);
    }
}
