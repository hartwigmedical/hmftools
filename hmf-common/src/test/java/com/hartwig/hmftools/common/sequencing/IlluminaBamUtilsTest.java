package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.getReadNameAttributes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.IlluminaReadNameAttributes;

import org.junit.Test;

public class IlluminaBamUtilsTest
{
    @Test
    public void testGetReadNameAttributes()
    {
        String readName = "A00624:8:HHKYHDSXX:4:2140:30409:19429";
        IlluminaReadNameAttributes expectedAttributes = new IlluminaReadNameAttributes("A00624", 8, "HHKYHDSXX", 4, 2140, 30409, 19429);
        IlluminaReadNameAttributes actualAttributes = getReadNameAttributes(readName);

        assertEquals(expectedAttributes, actualAttributes);
    }

    @Test
    public void testGetReadNameAttributesInvalidReadName()
    {
        String readName = "READ_001";
        IlluminaReadNameAttributes readNameAttributes = getReadNameAttributes(readName);

        assertNull(readNameAttributes);
    }

    @Test
    public void testGetReadNameAttributesInvalidRunId()
    {
        String readName = "A00624:X:HHKYHDSXX:4:2140:30409:19429";
        IlluminaReadNameAttributes readNameAttributes = getReadNameAttributes(readName);

        assertNull(readNameAttributes);
    }
}
