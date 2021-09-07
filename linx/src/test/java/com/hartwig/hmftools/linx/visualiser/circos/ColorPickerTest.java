package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.awt.Color;

import org.junit.Test;

public class ColorPickerTest
{
    @Test
    public void testContigColours() {
        Color chr1Color = ColorPicker.contigColour("1");
        assertNotEquals(Color.black, chr1Color);

        Color chr1ColorV38 = ColorPicker.contigColour("chr1");
        assertEquals(chr1Color, chr1ColorV38);

        Color unknownColor = ColorPicker.contigColour("MT");
        assertEquals(Color.black, unknownColor);
    }
}
