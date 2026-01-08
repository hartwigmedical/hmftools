package com.hartwig.hmftools.compar.image;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.ThresholdData;
import com.hartwig.hmftools.compar.common.ThresholdType;

public abstract class AbstractImageData implements ComparableItem
{
    public final String Name;
    public final String Path;
    public final BufferedImage Image;

    public static final String FLD_DIMENSION_MISMATCH = "DimensionMismatch";
    public static final String FLD_PIXEL_DIFF = "PixelDiff";

    public static final ThresholdData DEFAULT_IMAGE_THRESHOLD = new ThresholdData(
            ThresholdType.PERCENT, Double.NaN, 0);

    public AbstractImageData(String name, String path)
    {
        Name = name;
        Path = path;
        Image = loadImage(path);
    }

    public String getBasename(){ return new File(Path).getName(); }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean isPass() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final AbstractImageData otherImageData = (AbstractImageData) other;
        return Name.equals(otherImageData.Name);
    }

    @Override
    public String key() { return getBasename(); }

    @Override
    public List<String> displayValues() { return List.of(); }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final AbstractImageData otherImageData = (AbstractImageData) other;
        BufferedImage otherImage = otherImageData.Image;

        final List<String> diffs = Lists.newArrayList();

        if(Image.getWidth() != otherImage.getWidth() || Image.getHeight() != otherImage.getHeight())
        {
            String diffString = format("%s(%dx%d/%dx%d)",
                    FLD_DIMENSION_MISMATCH,
                    Image.getWidth(), Image.getHeight(),
                    otherImage.getWidth(), otherImage.getHeight()
            );
            diffs.add(diffString);
        }
        else
        {
            int totalPixels = Image.getWidth() * Image.getHeight();

            int absDiff = countDifferingPixels(Image, otherImage);
            double relDiff = (double) absDiff / totalPixels;

            ThresholdData threshold = thresholds.isFieldRegistered(Name)
                    ? thresholds.getThreshold(Name)
                    : DEFAULT_IMAGE_THRESHOLD;

            if(hasDiff(absDiff, relDiff, threshold))
            {
                String diffString = format("%s(%.3f=%d/%d)", FLD_PIXEL_DIFF, relDiff, absDiff, totalPixels);
                diffs.add(diffString);
            }
        }

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    public BufferedImage loadImage(String path)
    {
        try
        {
            BufferedImage image = ImageIO.read(new File(path));

            if(image == null)
            {
                CMP_LOGGER.error("failed to load image with unsupported format: {}", path);
                return null;
            }

            return image;
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load image: {}", path, e);
            return null;
        }
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

    public boolean hasDiff(int absDiff, double relDiff, ThresholdData threshold)
    {
        boolean hasAbsDiff = absDiff > threshold.AbsoluteDiff;
        boolean hasRelDiff = relDiff > threshold.PercentDiff;

        if(threshold.Type == ThresholdType.ABSOLUTE_AND_PERCENT)
        {
            return hasAbsDiff && hasRelDiff;
        }

        if(threshold.Type == ThresholdType.ABSOLUTE)
        {
            return hasAbsDiff;
        }

        return hasRelDiff;
    }
}
