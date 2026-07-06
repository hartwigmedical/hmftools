package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.awt.image.BufferedImage;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class ImageField implements Field
{
    public final String name;
    public final Function<ComparableItem, BufferedImage> extractValue;
    public final boolean isCompared;
    public final Double absoluteThreshold;
    public final Double percentThreshold;

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
    public String displayValue(final ComparableItem item)
    {
        return "";
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        BufferedImage oldImage = extractValue.apply(oldItem);
        BufferedImage newImage = extractValue.apply(newItem);

        int absDiff = countDifferingPixels(oldImage, newImage);
        if(absDiff == 0)
        {
            return false;
        }

        double relDiff = (double) absDiff / max(countTotalPixels(oldImage), countTotalPixels(newImage));

        boolean satisfiesAbsDiff = absoluteThreshold == null || absDiff > absoluteThreshold;
        boolean satisfiesRelDiff = percentThreshold == null || relDiff > percentThreshold;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }

    public List<String> determineDiffs(ComparableItem oldItem, ComparableItem newItem)
    {
        if(hasDiff(oldItem, newItem))
        {
            BufferedImage oldImage = extractValue.apply(oldItem);
            BufferedImage newImage = extractValue.apply(newItem);
            int absDiff = countDifferingPixels(oldImage, newImage);
            int totalPixels = max(countTotalPixels(oldImage), countTotalPixels(newImage));
            double relDiff = (double) absDiff / totalPixels;
            return List.of(format("%s(%.3f=%d/%d)", name(), relDiff, absDiff, totalPixels));
        }
        else
        {
            return Collections.emptyList();
        }
    }

    private static int countTotalPixels(BufferedImage image)
    {
        return image.getWidth() * image.getHeight();
    }

    private static int countDifferingPixels(BufferedImage image1, BufferedImage image2)
    {
        int diffCount = 0;

        for(int y = 0; y < image1.getHeight(); y++)
        {
            for(int x = 0; x < image1.getWidth(); x++)
            {
                if(image1.getRGB(x, y) != image2.getRGB(x, y))
                {
                    diffCount++;
                }
            }
        }

        return diffCount;
    }
}
