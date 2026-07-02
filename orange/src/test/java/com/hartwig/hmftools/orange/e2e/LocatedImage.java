package com.hartwig.hmftools.orange.e2e;

import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;

public record LocatedImage(Rectangle2D bounds, BufferedImage image)
{
    public double bottomYInPdfCoordinate()
    {
        return bounds.getY();
    }

    public double topYInPdfCoordinate()
    {
        return bounds.getY() + bounds.getHeight();
    }
}
