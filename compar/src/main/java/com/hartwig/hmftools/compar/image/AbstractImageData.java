package com.hartwig.hmftools.compar.image;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.imageio.ImageIO;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.ThresholdData;
import com.hartwig.hmftools.compar.common.ThresholdType;

public abstract class AbstractImageData implements ComparableItem
{
    private final Map<String, ImageRecord> Images;

    public static final String FLD_DIMENSION_MISMATCH = "DimensionMismatch";
    public static final String FLD_PIXEL_DIFF = "PixelDiff";

    public static final ThresholdData DEFAULT_PIXEL_DIFF_PERCENT_THRESHOLD = new ThresholdData(
            ThresholdType.PERCENT, Double.NaN, 0);

    public AbstractImageData(Map<String, String> imageNamePathMap)
    {
        LinkedHashMap<String, ImageRecord> images = Maps.newLinkedHashMap();

        for(String imageName : imageNamePathMap.values())
        {
            String imagePath = imageNamePathMap.get(imageName);
            BufferedImage image = loadImage(imagePath);

            if(image != null)
            {
                ImageRecord imageRecord = new ImageRecord(imageName, imagePath, image);
                images.put(imageName, imageRecord);
            }
        }

        Images = images;
    }

    public AbstractImageData(String imageName, String imagePath)
    {
        this(Map.of(imageName, imagePath));
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean isPass() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();

        for(String imageName : Images.keySet())
        {
            BufferedImage image = Images.get(imageName).Image;
            values.add(format("%s %dx%d", imageName, image.getWidth(), image.getHeight()));
        }

        return values;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final AbstractImageData otherData = (AbstractImageData) other;
        final List<String> diffs = Lists.newArrayList();

        for(String imageName : Images.keySet())
        {
            BufferedImage image = Images.get(imageName).Image;
            BufferedImage otherImage = otherData.Images.get(imageName).Image;

            String imageBasenames = Stream.of(Images.get(imageName).getBasename(), otherData.Images.get(imageName).getBasename())
                    .distinct()
                    .collect(Collectors.joining(", "));
            String imageMetadata = imageName + ": " + imageBasenames;

            if(image.getWidth() != otherImage.getWidth() || image.getHeight() != otherImage.getHeight())
            {
                String diffString = format("Image(%s) %s(%dx%d/%dx%d)",
                        imageMetadata, FLD_DIMENSION_MISMATCH,
                        image.getWidth(), image.getHeight(),
                        otherImage.getWidth(), otherImage.getHeight()
                );

                diffs.add(diffString);
            }
            else
            {
                int totalPixels = image.getWidth() * image.getHeight();

                int absDiff = countDifferingPixels(image, otherImage);
                double relDiff = (double) absDiff / totalPixels;

                ThresholdData threshold = thresholds.isFieldRegistered(imageName)
                        ? thresholds.getThreshold(imageName)
                        : DEFAULT_PIXEL_DIFF_PERCENT_THRESHOLD;

                if(hasDiff(absDiff, relDiff, threshold))
                {
                    String diffString = format("Image(%s) %s(%.3f=%d/%d)",
                            imageMetadata,
                            FLD_PIXEL_DIFF, relDiff, absDiff, totalPixels
                    );

                    diffs.add(diffString);
                }
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
