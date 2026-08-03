package com.hartwig.hmftools.compar.common;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.awt.Color;
import java.awt.image.BufferedImage;

import com.hartwig.hmftools.compar.common.field.BoolFieldValue;
import com.hartwig.hmftools.compar.common.field.DoubleFieldValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.IntFieldValue;
import com.hartwig.hmftools.compar.common.field.PixelFieldValue;
import com.hartwig.hmftools.compar.common.field.StringFieldValue;
import com.hartwig.hmftools.compar.common.field.ThresholdFieldCheck;

import org.junit.Test;

public class FieldValueTest
{
    private static final String FIELD_NAME = "Field";

    @Test
    public void testStringFieldValue()
    {
        StringFieldValue fieldValue1 = new StringFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), "value1");

        StringFieldValue fieldValue2 = new StringFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), "value1");

        assertFalse(fieldValue1.hasDifference(fieldValue2));

        fieldValue2 = new StringFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), "value2");

        assertTrue(fieldValue1.hasDifference(fieldValue2));
    }

    @Test
    public void testBoolFieldValue()
    {
        BoolFieldValue fieldValue1 = new BoolFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), true);

        BoolFieldValue fieldValue2 = new BoolFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), true);

        assertFalse(fieldValue1.hasDifference(fieldValue2));

        fieldValue2 = new BoolFieldValue(
                new FieldInfo(FIELD_NAME, new FieldCheck(true), null), false);

        assertTrue(fieldValue1.hasDifference(fieldValue2));
    }

    @Test
    public void testIntFieldValue()
    {
        ThresholdFieldCheck thresholdFieldCheck = new ThresholdFieldCheck(true, 10.0, 0.2);

        IntFieldValue fieldValue1 = new IntFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 80);

        IntFieldValue fieldValue2 = new IntFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 99);

        assertFalse(fieldValue1.hasDifference(fieldValue2));

        fieldValue2 = new IntFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 101);

        assertTrue(fieldValue1.hasDifference(fieldValue2));
    }

    @Test
    public void testDoubleFieldValue()
    {
        ThresholdFieldCheck thresholdFieldCheck = new ThresholdFieldCheck(true, 10.0, 0.2);

        DoubleFieldValue fieldValue1 = new DoubleFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 80);

        DoubleFieldValue fieldValue2 = new DoubleFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 99);

        assertFalse(fieldValue1.hasDifference(fieldValue2));

        fieldValue2 = new DoubleFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), 101);

        assertTrue(fieldValue1.hasDifference(fieldValue2));
    }

    private static final int IMAGE_WIDTH = 10;
    private static final int IMAGE_HEIGHT = 10;
    private static final int TOTAL_PIXELS = IMAGE_WIDTH * IMAGE_HEIGHT;

    // creates a white image with the first differingPixelCount pixels set to black
    private static BufferedImage createImage(final int differingPixelCount)
    {
        BufferedImage image = new BufferedImage(IMAGE_WIDTH, IMAGE_HEIGHT, BufferedImage.TYPE_INT_RGB);
        int white = Color.WHITE.getRGB();
        int black = Color.BLACK.getRGB();

        int pixelIndex = 0;
        for(int y = 0; y < IMAGE_HEIGHT; y++)
        {
            for(int x = 0; x < IMAGE_WIDTH; x++)
            {
                image.setRGB(x, y, pixelIndex < differingPixelCount ? black : white);
                pixelIndex++;
            }
        }
        return image;
    }

    @Test
    public void testPixelFieldValue()
    {
        ThresholdFieldCheck thresholdFieldCheck = new ThresholdFieldCheck(true, 1.0, 0.1);

        PixelFieldValue fieldValue1 = new PixelFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), createImage(10));

        PixelFieldValue fieldValue2 = new PixelFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), createImage(10));

        assertFalse(fieldValue1.hasDifference(fieldValue2));

        fieldValue2 = new PixelFieldValue(
                new FieldInfo(FIELD_NAME, thresholdFieldCheck, null), createImage(30));

        assertTrue(fieldValue1.hasDifference(fieldValue2));
    }

}
