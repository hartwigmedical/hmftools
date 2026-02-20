package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.common.ThresholdData;
import com.hartwig.hmftools.compar.common.ThresholdType;

public abstract class ComparableImage implements ComparableItem
{
    public final String Name;
    public final String Path;
    public final BufferedImage Image;

    public static final String FLD_DIMENSIONS = "Dimensions";
    public static final String FLD_PIXELS = "Pixels";

    public static final ThresholdData DEFAULT_IMAGE_THRESHOLD = new ThresholdData(ThresholdType.PERCENT, Double.NaN, 0);

    public ComparableImage(String name, String path)
    {
        Name = name;
        Path = path;
        Image = loadImage(path);
    }

    @VisibleForTesting
    protected ComparableImage(String name, BufferedImage image)
    {
        Name = name;
        Path = null;
        Image = image;
    }

    public String getBasename(){ return new File(Path).getName(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final ComparableImage otherImageData = (ComparableImage) other;
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
        final ComparableImage otherImageData = (ComparableImage) other;
        BufferedImage otherImage = otherImageData.Image;

        // Handle image load failures
        if(Image == null && otherImage == null)
            return new Mismatch(this, other, MismatchType.INVALID_BOTH, List.of());

        if(Image == null)
            return new Mismatch(this, other, MismatchType.INVALID_REF, List.of());

        if(otherImage == null)
            return new Mismatch(this, other, MismatchType.INVALID_NEW, List.of());

        // Check dimensions
        if(Image.getWidth() != otherImage.getWidth() || Image.getHeight() != otherImage.getHeight())
        {
            String diffString = formDimensionMismatchString(Image.getWidth(), Image.getHeight(), otherImage.getWidth(), otherImage.getHeight());
            return new Mismatch(this, other, MismatchType.VALUE, List.of(diffString));
        }

        // Check pixels
        int totalPixels = Image.getWidth() * Image.getHeight();
        int absDiff = countDifferingPixels(Image, otherImage);
        double relDiff = (double) absDiff / totalPixels;
        ThresholdData threshold = thresholds.isFieldRegistered(Name) ? thresholds.getThreshold(Name) : DEFAULT_IMAGE_THRESHOLD;

        if(hasDiff(absDiff, relDiff, threshold))
        {
            String diffString = formPixelMismatchString(relDiff, absDiff, totalPixels);
            return new Mismatch(this, other, MismatchType.VALUE, List.of(diffString));
        }

        if(includeMatches)
            return new Mismatch(this, other, MismatchType.FULL_MATCH, List.of());

        return null;
    }

    @VisibleForTesting
    protected static String formDimensionMismatchString(int width, int height, int otherWidth, int otherHeight)
    {
        return String.format("%s(%dx%d/%dx%d)", FLD_DIMENSIONS, width, height, otherWidth, otherHeight);
    }

    @VisibleForTesting
    protected static String formPixelMismatchString(double relDiff, int absDiff, int totalPixels)
    {
        return String.format("%s(%.3f=%d/%d)", FLD_PIXELS, relDiff, absDiff, totalPixels);
    }

    public BufferedImage loadImage(String path)
    {
        try
        {
            BufferedImage image = ImageIO.read(new File(path));

            if(image == null)
            {
                CMP_LOGGER.warn("failed to load image with unsupported format: {}", path);
                return null;
            }

            return image;
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("failed to load image: {}", path, e);
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

    private boolean hasDiff(int absDiff, double relDiff, ThresholdData threshold)
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
