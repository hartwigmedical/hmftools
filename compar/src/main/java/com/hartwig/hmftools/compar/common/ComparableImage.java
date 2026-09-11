package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.ImageComparer.FLD_DIMENSIONS;
import static com.hartwig.hmftools.compar.common.ImageComparer.FLD_PIXELS;
import static com.hartwig.hmftools.compar.common.field.FieldInfo.findField;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.PixelFieldValue;

public abstract class ComparableImage extends ComparableItem
{
    public final String Name;
    public final String Path;
    public final BufferedImage Image;

    public ComparableImage(final String name, final String path, final List<FieldInfo> fields)
    {
        Name = name;
        Path = path;
        Image = loadImage(path);

        addStringValue(FLD_DIMENSIONS, dimensionString(), fields);
        mValues.put(FLD_PIXELS, new PixelFieldValue(findField(FLD_PIXELS, fields), Image));
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
