package com.hartwig.hmftools.esvee.prep;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testReadTrimmer()
    {
        ReadIdTrimmer trimmer = new ReadIdTrimmer(false);

        String readId = "READ_01";
        assertEquals(readId, trimmer.trim(readId));

        trimmer = new ReadIdTrimmer(true);

        // illumina style: A00260:251:HLYGFDSXY:1:1673:32280:4946
        readId = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
        assertEquals("1:1673:32280:4946", trimmer.trim(readId));

        readId = "A00260:251:HLYGFDSXY:1:2:3:4";
        assertEquals("1:2:3:4", trimmer.trim(readId));

        // disable on a read with different delims or too short
        readId = "A0:25:HL:1:2:3:4";
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());

        trimmer = new ReadIdTrimmer(true); // reset

        readId = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
        assertEquals("1:1673:32280:4946", trimmer.trim(readId));

        readId = "A00260:251:HLYGFDSX:1:2:3:4"; // diff delim position
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());

        // other formats are not currently supported
        trimmer = new ReadIdTrimmer(true); // reset

        readId = "011852_2-UGAv3-2-1333458495";
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());
    }
}
