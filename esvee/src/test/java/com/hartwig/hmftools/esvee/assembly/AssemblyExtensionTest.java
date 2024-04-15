package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RefBaseAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

import htsjdk.samtools.util.StringUtil;

public class AssemblyExtensionTest
{
    @Test
    public void testRefBaseExtension()
    {
        // index           0         10        20        30        40         50
        // pos             30        40        50        60        70         80
        //                 01234567890012345678901234567890123456789012345678901234567890
        String refBases = "AAACCCGGGTTTAACCGGTTAAAAACCCCCGGGGGTTTTTAACCGGTTAAACCCGGGTTTAA";

        Junction posJunction = new Junction(CHR_1, 60, POS_ORIENT);

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(refBases.length());

        String existingRefBases = refBases.substring(20, 41);
        JunctionAssembly assembly = new JunctionAssembly(posJunction, existingRefBases.getBytes(), baseQuals, 10);

        RefBaseAssembly refBaseAssembly = new RefBaseAssembly(assembly, 40, null);

        assertEquals(21, refBaseAssembly.baseLength());
        String refAssemblyBases = StringUtil.bytesToString(refBaseAssembly.bases(), 10, 11);
        assertEquals("AAAAACCCCCG", refAssemblyBases);

        // now test adding reads to this section
        int permittedMismatches = 1;
        int requiredOverlap = 5;

        // outside the bounds of the extended region
        Read read1 = createRead("READ_01", 39, refBases.substring(9, 29), "20M");
        assertFalse(refBaseAssembly.checkAddRead(read1, JUNCTION_MATE, 1, 10));

        Read read2 = createRead("READ_02", 40, refBases.substring(10, 25), "15M");
        assertTrue(refBaseAssembly.checkAddRead(read2, JUNCTION_MATE, permittedMismatches, requiredOverlap));

        // with 1 mismatch
        Read read3 = createRead(
                "READ_03", 45, refBases.substring(15, 20) + "G" + refBases.substring(21, 30), "15M");

        assertTrue(refBaseAssembly.checkAddRead(read3, DISCORDANT,permittedMismatches, requiredOverlap));
        assertEquals(1, refBaseAssembly.readMismatchesCount());

        // with 2 mismatches
        Read read4 = createRead(
                "READ_04", 45, refBases.substring(15, 19) + "GG" + refBases.substring(21, 30), "15M");

        assertFalse(refBaseAssembly.checkAddRead(read4, DISCORDANT, permittedMismatches, requiredOverlap));

        Junction negJunction = new Junction(CHR_1, 60, NEG_ORIENT);

        assembly = new JunctionAssembly(negJunction, existingRefBases.getBytes(), baseQuals, 10); // was 50 -> 70

        refBaseAssembly = new RefBaseAssembly(assembly, 80, null);

        assertEquals(21, refBaseAssembly.baseLength());
        refAssemblyBases = StringUtil.bytesToString(refBaseAssembly.bases(), 0, 11);
        assertEquals("GGGGGTTTTTA", refAssemblyBases);

        // repeat tests for reads on the ref side of a negative junction
        read1 = createRead("READ_01", 81, refBases.substring(51, 56), "5M");
        assertFalse(refBaseAssembly.checkAddRead(read1, DISCORDANT,permittedMismatches, requiredOverlap));

        read2 = createRead("READ_02", 65, refBases.substring(35, 50), "15M");
        assertTrue(refBaseAssembly.checkAddRead(read2, DISCORDANT,permittedMismatches, requiredOverlap));

        // with 1 mismatch
        read3 = createRead(
                "READ_03", 70, refBases.substring(40, 45) + "C" + refBases.substring(46, 50), "10M");

        assertTrue(refBaseAssembly.checkAddRead(read3, DISCORDANT, permittedMismatches, requiredOverlap));
        assertEquals(1, refBaseAssembly.readMismatchesCount());

        // with 2 mismatches
        read4 = createRead(
                "READ_04", 70, refBases.substring(40, 44) + "CC" + refBases.substring(46, 50), "10M");

        assertFalse(refBaseAssembly.checkAddRead(read4, DISCORDANT, permittedMismatches, requiredOverlap));
    }
}
