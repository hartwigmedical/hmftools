package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.awt.image.BufferedImage;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class ImageField implements Field
{
    private final String name;
    private final Function<ComparableItem, BufferedImage> extractValue;
    private final boolean isCompared;
    private final Double absoluteThreshold;
    private final Double percentThreshold;

    public ImageField(final String name, final Function<ComparableItem, BufferedImage> extractValue, final boolean isCompared,
            final Double absoluteThreshold, final Double percentThreshold)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.isCompared = isCompared;
        this.absoluteThreshold = absoluteThreshold;
        this.percentThreshold = percentThreshold;
    }

    @Override
    public String name()
    {
        return name;
    }

    @Override
    public boolean isCompared()
    {
        return isCompared;
    }

    @Override
    public Double absoluteThreshold()
    {
        return absoluteThreshold;
    }

    @Override
    public Double percentThreshold()
    {
        return percentThreshold;
    }

    @Override
    public String type()
    {
        return "image";
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return "";
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return !determineDiffs(oldItem, newItem).isEmpty();
    }

    @Override
    public List<String> determineDiffs(ComparableItem oldItem, ComparableItem newItem)
    {
        BufferedImage oldImage = extractValue.apply(oldItem);
        BufferedImage newImage = extractValue.apply(newItem);

        int absDiff = countDifferingPixels(oldImage, newImage);
        if(absDiff == 0)
        {
            return Collections.emptyList();
        }

        int totalPixels = countSpanningPixels(oldImage, newImage);
        double relDiff = (double) absDiff / totalPixels;

        boolean satisfiesAbsDiff = absoluteThreshold == null || absDiff > absoluteThreshold;
        boolean satisfiesRelDiff = percentThreshold == null || relDiff > percentThreshold;

        boolean hasDiff = satisfiesAbsDiff && satisfiesRelDiff;

        if(hasDiff)
        {
            return List.of(format("%s(%.3f=%d/%d)", name(), relDiff, absDiff, totalPixels));
        }
        else
        {
            return Collections.emptyList();
        }
    }

    private static int countSpanningPixels(BufferedImage oldImage, BufferedImage newImage)
    {
        return max(oldImage.getHeight(), newImage.getHeight()) * max(oldImage.getWidth(), newImage.getWidth());
    }

    private static int countDifferingPixels(BufferedImage oldImage, BufferedImage newImage)
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
