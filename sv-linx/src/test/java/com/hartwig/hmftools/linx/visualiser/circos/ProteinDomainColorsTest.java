package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ProteinDomainColorsTest
{
    @Test
    public void testHue() {

        assertEquals(21 + 90, 360 * ProteinDomainColors.hue(0, 1), 0.1);

        assertEquals(21 + 90, 360 * ProteinDomainColors.hue(0, 2), 0.1);
        assertEquals(21 + 270, 360 * ProteinDomainColors.hue(1, 2), 0.1);

        assertEquals(21 + 60, 360 * ProteinDomainColors.hue(0, 3), 0.1);
        assertEquals(21 + 180 + 60, 360 * ProteinDomainColors.hue(1, 3), 0.1);
        assertEquals(21 + 120, 360 * ProteinDomainColors.hue(2, 3), 0.1);

        assertEquals(21 + 60, 360 * ProteinDomainColors.hue(0, 4), 0.1);
        assertEquals(21 + 180 + 60, 360 * ProteinDomainColors.hue(1, 4), 0.1);
        assertEquals(21 + 120, 360 * ProteinDomainColors.hue(2, 4), 0.1);
        assertEquals(21 + 180 + 120, 360 * ProteinDomainColors.hue(3, 4), 0.1);

    }

}
