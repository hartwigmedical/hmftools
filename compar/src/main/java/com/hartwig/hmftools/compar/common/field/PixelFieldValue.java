package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.awt.image.BufferedImage;
import java.util.List;

import com.google.common.collect.Lists;

public class PixelFieldValue extends FieldValue
{
    public final BufferedImage Image;

    public PixelFieldValue(final FieldInfo field, final BufferedImage image)
    {
        super(field);
        Image = image;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        List<String> diffs = Lists.newArrayList();
        addDifferences(other, diffs);
        return !diffs.isEmpty();
    }

    private void addDifferences(final FieldValue other, final List<String> diffs)
    {
        PixelFieldValue otherValue = (PixelFieldValue)other;

        BufferedImage otherImage = otherValue.Image;

        int absDiff = countDifferingPixels(Image, otherImage);
        if(absDiff == 0)
            return;

        int totalPixels = countSpanningPixels(otherImage, Image);
        double percDiff = (double) absDiff / totalPixels;

        ThresholdFieldCheck thresholdFieldCheck = (ThresholdFieldCheck)Field.FieldCheck;

        boolean satisfiesAbsDiff = thresholdFieldCheck.AbsoluteDiff == null || absDiff > thresholdFieldCheck.AbsoluteDiff;
        boolean satisfiesRelDiff = thresholdFieldCheck.PercentageDiff == null || percDiff > thresholdFieldCheck.PercentageDiff;

        boolean hasDiff = satisfiesAbsDiff && satisfiesRelDiff;

        if(hasDiff)
        {
            diffs.add(format("%s(%.3f=%d/%d)", name(), percDiff, absDiff, totalPixels));
        }
    }

    @Override
    public String toString() { return format("%s=%s", name(), formatString()); }

    @Override
    public String displayValue() { return ""; }

    private static int countSpanningPixels(final BufferedImage oldImage, final BufferedImage newImage)
    {
        return max(oldImage.getHeight(), newImage.getHeight()) * max(oldImage.getWidth(), newImage.getWidth());
    }

    private static int countDifferingPixels(final BufferedImage oldImage, final BufferedImage newImage)
    {
        int diffCount = 0;

        int minHeight = min(oldImage.getHeight(), newImage.getHeight());
        int maxHeight = max(oldImage.getHeight(), newImage.getHeight());
        int minWidth = min(oldImage.getWidth(), newImage.getWidth());
        int maxWidth = max(oldImage.getWidth(), newImage.getWidth());

        for(int y = 0; y < minHeight; y++)
        {
            for(int x = 0; x < minWidth; x++)
            {
                if(oldImage.getRGB(x, y) != newImage.getRGB(x, y))
                {
                    diffCount++;
                }
            }
        }

        // count pixels in one image but not the other as differences
        diffCount += maxHeight * maxWidth - minHeight * minWidth;

        return diffCount;
    }
}
