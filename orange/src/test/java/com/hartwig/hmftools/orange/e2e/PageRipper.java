package com.hartwig.hmftools.orange.e2e;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.List;

import com.google.common.base.Preconditions;

import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.text.PDFTextStripperByArea;
import org.jetbrains.annotations.NotNull;
import org.junit.Assert;

public class PageRipper
{
    private final PDPage page;
    private final int pageHeight;
    private final int pageWidth;
    private final List<LocatedImage> images;

    public PageRipper(final PDPage page) throws IOException
    {
        this.page = page;
        pageHeight = (int) page.getMediaBox().getHeight();
        pageWidth = (int) page.getMediaBox().getWidth();
        PDFImageLocator locator = new PDFImageLocator(page);
        images = locator.locateImages();
    }

    public String[] getLinesInRectangle(PositionInPage topLeft, PositionInPage bottomRight) throws IOException
    {
        Preconditions.checkArgument(topLeft.fractionAcrossPage < bottomRight.fractionAcrossPage);
        Preconditions.checkArgument(topLeft.fractionDownPage < bottomRight.fractionDownPage);
        int topLeftX = (int) (topLeft.fractionAcrossPage * pageWidth);
        int topLeftY = (int) (topLeft.fractionDownPage * pageHeight);
        int width = (int) ((bottomRight.fractionAcrossPage - topLeft.fractionAcrossPage) * pageWidth);
        int height = (int) ((bottomRight.fractionDownPage - topLeft.fractionDownPage) * pageHeight);
        Rectangle region = new Rectangle(topLeftX, topLeftY, width, height);

        PDFTextStripperByArea byArea = new PDFTextStripperByArea();
        byArea.setSortByPosition(true);
        byArea.addRegion("region", region);
        byArea.extractRegions(page);
        return byArea.getTextForRegion("region").split("\n");
    }

    public int numberOfImages()
    {
        return images.size();
    }

    public LocatedImage checkHasImageWithinBoundsOfGivenSize(RectangleInPage expectedBounds, int expectedWidth, int expectedHeight)
    {
        LocatedImage imageWithinBounds = getUniqueImageWithinBounds(expectedBounds);

        Assert.assertEquals(expectedWidth, imageWithinBounds.image().getWidth());
        Assert.assertEquals(expectedHeight, imageWithinBounds.image().getHeight());

        return imageWithinBounds;
    }

    public void checkHasImageWithinBoundsOfGivenSizeAndColor(RectangleInPage expectedBounds,
            int expectedWidth, int expectedHeight, Color expectedColor)
    {
        LocatedImage locatedImage = checkHasImageWithinBoundsOfGivenSize(expectedBounds, expectedWidth, expectedHeight);
        BufferedImage imgData = locatedImage.image();
        int referenceRgb = imgData.getRGB(0, 0);

        boolean isSingleColor = true;
        outerLoop:
        for(int y = 0; y < imgData.getHeight(); y++)
        {
            for(int x = 0; x < imgData.getWidth(); x++)
            {
                if(imgData.getRGB(x, y) != referenceRgb)
                {
                    isSingleColor = false;
                    break outerLoop;
                }
            }
        }

        Assert.assertTrue("Image within bounds is not single color", isSingleColor);
        Color actualColor = new Color(referenceRgb, true);
        Assert.assertEquals(expectedColor.getRed(), actualColor.getRed());
        Assert.assertEquals(expectedColor.getGreen(), actualColor.getGreen());
        Assert.assertEquals(expectedColor.getBlue(), actualColor.getBlue());
    }

    @NotNull
    private LocatedImage getUniqueImageWithinBounds(final RectangleInPage expectedBounds)
    {
        LocatedImage imageWithinBounds = null;
        for(LocatedImage locatedImage : images)
        {
            RectangleInPage actualBounds = convertRectangle(locatedImage);
            if(expectedBounds.contains(actualBounds))
            {
                if(imageWithinBounds != null)
                {
                    throw new AssertionError("Multiple images found within bounds " + expectedBounds);
                }
                imageWithinBounds = locatedImage;
            }
        }
        Assert.assertNotNull(imageWithinBounds);
        return imageWithinBounds;
    }

    private RectangleInPage convertRectangle(LocatedImage locatedImage)
    {
        double topLeftX = locatedImage.bounds().getX() / pageWidth;
        double topLeftY = (pageHeight - locatedImage.topYInPdfCoordinate()) / pageHeight;
        double fractionalWidth = locatedImage.bounds().getWidth() / pageWidth;
        double fractionalHeight = locatedImage.bounds().getHeight() / pageHeight;
        double bottomRightX = topLeftX + fractionalWidth;
        double bottomRightY = topLeftY + fractionalHeight;

        return new RectangleInPage(topLeftX, topLeftY, bottomRightX, bottomRightY);
    }
}
