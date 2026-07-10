package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public abstract class ComparableImage implements ComparableItem
{
    public final String Name;
    public final String Path;
    public final BufferedImage Image;

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

    public String dimensionString()
    {
        return String.format("%dx%d", Image.getWidth(), Image.getHeight());
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
}
