package com.hartwig.hmftools.compar;

import static org.junit.Assert.assertEquals;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Test;

public class ComparableImageTest
{
    private static final int TEST_IMAGE_WIDTH = 10;
    private static final int TEST_IMAGE_HEIGHT = 10;
    private static final DiffThresholds TEST_THRESHOLDS = new DiffThresholds();

    private static class TestImageData extends ComparableImage
    {
        public TestImageData(String name, BufferedImage image) { super(name, image); }

        @Override
        public Category category() { return Category.CUPPA_IMAGE; }
    }

    private TestImageData createTestBlackAndWhiteImage(int width, int height, double blackProportion)
    {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        int blackWhiteHeightMidpoint = (int) (height * blackProportion);

        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width, blackWhiteHeightMidpoint);

        g2d.setColor(Color.BLACK);
        g2d.fillRect(0, blackWhiteHeightMidpoint, width, height-blackWhiteHeightMidpoint);

        g2d.dispose();
        return new TestImageData("test", image);
    }

    private TestImageData createTestBlackAndWhiteImage(double blackProportion)
    {
        return createTestBlackAndWhiteImage(TEST_IMAGE_WIDTH, TEST_IMAGE_HEIGHT, blackProportion);
    }

    @Test
    public void differentDimensionsProduceDimensionMismatch()
    {
        TestImageData whiteImage1 = createTestBlackAndWhiteImage(1, 1, 0);
        TestImageData whiteImage2 = createTestBlackAndWhiteImage(2, 2, 0);

        Mismatch mismatch = whiteImage1.findMismatch(whiteImage2, MatchLevel.DETAILED, TEST_THRESHOLDS, true);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals("DimensionMismatch(1x1/2x2)", mismatch.DiffValues().get(0));
    }

    @Test
    public void identicalWhiteImagesMatch()
    {
        TestImageData whiteImage1 = createTestBlackAndWhiteImage(0);
        TestImageData whiteImage2 = createTestBlackAndWhiteImage(0);

        Mismatch mismatch = whiteImage1.findMismatch(whiteImage2, MatchLevel.DETAILED, TEST_THRESHOLDS, true);

        assertEquals(MismatchType.FULL_MATCH, mismatch.MismatchType());
    }

    @Test
    public void whiteVsBlackImageHas100PercentMismatch()
    {
        TestImageData whiteImage = createTestBlackAndWhiteImage(0);
        TestImageData blackImage = createTestBlackAndWhiteImage(1);

        Mismatch mismatch = whiteImage.findMismatch(blackImage, MatchLevel.DETAILED, TEST_THRESHOLDS, true);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals("PixelMismatch(1.000=100/100)", mismatch.DiffValues().get(0));
    }

    @Test
    public void whiteVsHalfBlackHalfWhiteImageHas50PercentMismatch()
    {
        TestImageData whiteImage = createTestBlackAndWhiteImage(0);
        TestImageData mixedImage = createTestBlackAndWhiteImage(0.5);

        Mismatch mismatch = whiteImage.findMismatch(mixedImage, MatchLevel.DETAILED, TEST_THRESHOLDS, true);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals("PixelMismatch(0.500=50/100)", mismatch.DiffValues().get(0));
    }
}