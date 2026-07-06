package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.List;

import org.junit.Test;

public class ImageFieldTest
{
    private static final String FIELD_NAME = "ImageField";
    private static final int IMAGE_WIDTH = 10;
    private static final int IMAGE_HEIGHT = 10;
    private static final int TOTAL_PIXELS = IMAGE_WIDTH * IMAGE_HEIGHT;

    private static ImageField field(final Double absoluteThreshold, final Double percentThreshold)
    {
        return new ImageField(FIELD_NAME, i -> ((TestFieldItem<BufferedImage>) i).Value, true, absoluteThreshold, percentThreshold);
    }

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

    private static boolean hasDiff(final ImageField field, final int oldDifferingPixels, final int newDifferingPixels)
    {
        return field.hasDiff(new TestFieldItem<>(createImage(oldDifferingPixels)), new TestFieldItem<>(createImage(newDifferingPixels)));
    }

    @Test
    public void nameIsComparedAndThresholdsReflectConstructorArgs()
    {
        ImageField field = field(5., 0.2);
        assertEquals(FIELD_NAME, field.name());
        assertTrue(field.isCompared());
        assertEquals(5., field.absoluteThreshold(), 0.0);
        assertEquals(0.2, field.percentThreshold(), 0.0);
    }

    @Test
    public void displayValueIsAlwaysEmpty()
    {
        ImageField field = field(0., null);
        assertEquals("", field.displayValue(new TestFieldItem<>(createImage(50))));
    }

    @Test
    public void hasDiffIsFalseForIdenticalImages()
    {
        ImageField field = field(0., null);
        assertFalse(hasDiff(field, 0, 0));
        assertFalse(hasDiff(field, 30, 30));
    }

    @Test
    public void withoutThresholdsIdenticalImagesDoNotDiffButAnyDifferenceDoes()
    {
        ImageField field = field(null, null);
        assertFalse(hasDiff(field, 0, 0));
        assertTrue(hasDiff(field, 0, 1));
    }

    @Test
    public void absoluteThresholdOnlyRequiresDifferingPixelCountToExceedThreshold()
    {
        ImageField field = field(5., null);
        assertFalse(hasDiff(field, 0, 5));
        assertTrue(hasDiff(field, 0, 6));
    }

    @Test
    public void percentThresholdOnlyRequiresProportionToExceedThreshold()
    {
        ImageField field = field(null, 0.2);
        assertFalse(hasDiff(field, 0, 20));
        assertTrue(hasDiff(field, 0, 21));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceeded()
    {
        // absolute threshold lower than percent threshold's pixel-count equivalent (50): allows isolating "absolute exceeded only"
        ImageField absDominant = field(20., 0.5);
        assertFalse(hasDiff(absDominant, 0, 25)); // 25 > 20 (abs) but 25% not > 50% (percent) -> no diff
        assertTrue(hasDiff(absDominant, 0, 60)); // 60 > 20 (abs) and 60% > 50% (percent) -> diff

        // percent threshold's pixel-count equivalent (20) lower than absolute threshold (60): allows isolating "percent exceeded only"
        ImageField percentDominant = field(60., 0.2);
        assertFalse(hasDiff(percentDominant, 0, 25)); // 25% > 20% (percent) but 25 not > 60 (abs) -> no diff
    }

    @Test
    public void determineDiffsFormatsProportionAndPixelCounts()
    {
        ImageField field = field(null, 0.);
        List<String> diffs = field.determineDiffs(new TestFieldItem<>(createImage(0)), new TestFieldItem<>(createImage(50)));
        assertEquals(List.of(String.format("%s(%.3f=%d/%d)", FIELD_NAME, 0.5, 50, TOTAL_PIXELS)), diffs);
    }

    @Test
    public void determineDiffsIsEmptyWhenNoDiff()
    {
        ImageField field = field(null, 0.);
        assertTrue(field.determineDiffs(new TestFieldItem<>(createImage(30)), new TestFieldItem<>(createImage(30))).isEmpty());
    }
}
