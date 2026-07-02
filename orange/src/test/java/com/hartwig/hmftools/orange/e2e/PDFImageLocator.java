package com.hartwig.hmftools.orange.e2e;

import org.apache.pdfbox.contentstream.PDFGraphicsStreamEngine;
import org.apache.pdfbox.contentstream.PDFStreamEngine;
import org.apache.pdfbox.contentstream.operator.DrawObject;
import org.apache.pdfbox.contentstream.operator.Operator;
import org.apache.pdfbox.contentstream.operator.OperatorName;
import org.apache.pdfbox.cos.COSBase;
import org.apache.pdfbox.cos.COSName;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.graphics.PDXObject;
import org.apache.pdfbox.pdmodel.graphics.image.PDImage;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;
import org.apache.pdfbox.util.Matrix;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class PDFImageLocator extends PDFGraphicsStreamEngine
{

    private final List<LocatedImage> imagesOnPage = new ArrayList<>();

    // Constructor maps directly to PDFGraphicsStreamEngine constraints
    public PDFImageLocator(PDPage page)
    {
        super(page);
    }

    public List<LocatedImage> locateImages() throws IOException
    {
        imagesOnPage.clear();
        processPage(getPage()); // This starts the stream processing
        return new ArrayList<>(imagesOnPage);
    }

    // --- Obligatory Blank Overrides Required by PDFGraphicsStreamEngine ---
    @Override
    public void appendRectangle(Point2D p0, Point2D p1, Point2D p2, Point2D p3) throws IOException
    {
    }

    @Override
    public void drawImage(PDImage pdImage) throws IOException
    {
        // PDFBox's layout rendering subsystem stores image resources as PDImageXObject instances
        if(pdImage instanceof PDImageXObject)
        {
            PDImageXObject imageXObject = (PDImageXObject) pdImage;

            // This grabs the fully calculated context transformation matrix up to this exact draw execution
            Matrix ctm = getGraphicsState().getCurrentTransformationMatrix();

            // Transform the 4 vertices of the standard coordinate matrix
            Point2D p0 = ctm.transformPoint(0, 0);
            Point2D p1 = ctm.transformPoint(1, 0);
            Point2D p2 = ctm.transformPoint(0, 1);
            Point2D p3 = ctm.transformPoint(1, 1);

            // Calculate exact bounding boxes regardless of layout directions or mirroring
            double minX = Math.min(Math.min(p0.getX(), p1.getX()), Math.min(p2.getX(), p3.getX()));
            double maxX = Math.max(Math.max(p0.getX(), p1.getX()), Math.max(p2.getX(), p3.getX()));
            double minY = Math.min(Math.min(p0.getY(), p1.getY()), Math.min(p2.getY(), p3.getY()));
            double maxY = Math.max(Math.max(p0.getY(), p1.getY()), Math.max(p2.getY(), p3.getY()));

            double absoluteWidth = maxX - minX;
            double absoluteHeight = maxY - minY;

            Rectangle2D absoluteBounds = new Rectangle2D.Double(minX, minY, absoluteWidth, absoluteHeight);

            // Extract the buffered bitmap pixel array for color testing
            imagesOnPage.add(new LocatedImage(absoluteBounds, imageXObject.getImage()));
        }
    }

    @Override
    public void clip(int windingRule) throws IOException
    {
    }

    @Override
    public void moveTo(float x, float y) throws IOException
    {
    }

    @Override
    public void lineTo(float x, float y) throws IOException
    {
    }

    @Override
    public void curveTo(float x1, float y1, float x2, float y2, float x3, float y3) throws IOException
    {
    }

    @Override
    public Point2D getCurrentPoint() throws IOException
    {
        return new Point2D.Float(0, 0);
    }

    @Override
    public void closePath() throws IOException
    {
    }

    @Override
    public void endPath() throws IOException
    {
    }

    @Override
    public void strokePath() throws IOException
    {
    }

    @Override
    public void fillPath(int windingRule) throws IOException
    {
    }

    @Override
    public void fillAndStrokePath(int windingRule) throws IOException
    {
    }

    @Override
    public void shadingFill(COSName shadingName) throws IOException
    {
    }
}