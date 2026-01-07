package com.hartwig.hmftools.compar.image;

import java.awt.image.BufferedImage;
import java.io.File;

public class ImageRecord
{
    public final String Name;
    public final String Path;
    public final BufferedImage Image;

    public ImageRecord(String name, String path, BufferedImage image)
    {
        Name = name;
        Path = path;
        Image = image;
    }

    public String getBasename() { return new File(Path).getName(); }
}

