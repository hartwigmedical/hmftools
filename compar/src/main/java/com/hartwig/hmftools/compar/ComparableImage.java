package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import com.google.common.annotations.VisibleForTesting;

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
