package com.hartwig.hmftools.orange.report;

import java.awt.Color;
import java.io.IOException;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.font.PDFont;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;

import be.quodlibet.boxable.BaseTable;

public class DocumentContext
{
    private final PDDocument mDocument;
    private PDPage mCurrentPage;
    private PDRectangle mPageSize;

    private float mCursorY;

    private final float mMarginTop;
    private final float mMarginBottom;
    private final float mMarginLeft;
    private final float mMarginRight;

    private PageEventHandler mPageEventHandler;

    public DocumentContext(final PDDocument document, final PDRectangle pageSize,
            float marginTop, float marginBottom, float marginLeft, float marginRight)
    {
        mDocument = document;
        mPageSize = pageSize;
        mMarginTop = marginTop;
        mMarginBottom = marginBottom;
        mMarginLeft = marginLeft;
        mMarginRight = marginRight;
        mCurrentPage = null;
        mCursorY = 0;
    }

    public void setPageEventHandler(final PageEventHandler handler)
    {
        mPageEventHandler = handler;
    }

    public PDDocument document()
    {
        return mDocument;
    }

    public PDPage currentPage()
    {
        return mCurrentPage;
    }

    public PDRectangle pageSize()
    {
        return mPageSize;
    }

    public void setPageSize(final PDRectangle pageSize)
    {
        mPageSize = pageSize;
    }

    public float cursorY()
    {
        return mCursorY;
    }

    public float contentWidth()
    {
        return mPageSize.getWidth() - mMarginLeft - mMarginRight;
    }

    public float contentStartY()
    {
        return mPageSize.getHeight() - mMarginTop;
    }

    public float contentEndY()
    {
        return mMarginBottom;
    }

    public float marginLeft()
    {
        return mMarginLeft;
    }

    /**
     * Create a new page with the current page size and reset cursor.
     */
    public void newPage() throws IOException
    {
        newPage(mPageSize);
    }

    /**
     * Create a new page with the specified page size and reset cursor.
     */
    public void newPage(final PDRectangle pageSize)
    {
        mPageSize = pageSize;
        mCurrentPage = new PDPage(pageSize);
        mDocument.addPage(mCurrentPage);
        mCursorY = contentStartY();

        if(mPageEventHandler != null)
        {
            mPageEventHandler.onPageStart(mCurrentPage, mDocument);
        }
    }

    /**
     * Ensure we have a current page; create one if needed.
     */
    private void ensurePage() throws IOException
    {
        if(mCurrentPage == null)
        {
            newPage();
        }
    }

    /**
     * Check if the remaining space on the current page can fit the given height.
     */
    private boolean hasSpaceFor(float height)
    {
        return mCurrentPage != null && (mCursorY - height) >= contentEndY();
    }

    /**
     * Draw a Boxable BaseTable at the current cursor position.
     * The table is drawn and the cursor is advanced by the table height.
     * Boxable handles multi-page overflow internally if configured.
     */
    public void addTable(final BaseTable table) throws IOException
    {
        ensurePage();
        // draw() returns the final Y position after rendering
        float finalY = table.draw();
        // Update current page in case the table spanned multiple pages
        if(table.getCurrentPage() != mCurrentPage)
        {
            mCurrentPage = table.getCurrentPage();
            if(mPageEventHandler != null)
            {
                mPageEventHandler.onPageStart(mCurrentPage, mDocument);
            }
        }
        mCursorY = finalY - 5; // small spacing after table
    }

    /**
     * Create a Boxable BaseTable ready for row/cell additions.
     * The table starts at the current cursor Y position.
     */
    public BaseTable createTable(float width, float[] columnWidths) throws IOException
    {
        ensurePage();
        float yStart = mCursorY;
        float yStartNewPage = contentStartY();
        float bottomMargin = contentEndY();

        BaseTable table = new BaseTable(yStart, yStartNewPage, bottomMargin, width, mMarginLeft, mDocument, mCurrentPage, true, true);
        return table;
    }

    /**
     * Create a simple Boxable BaseTable with full content width.
     */
    public BaseTable createTable() throws IOException
    {
        return createTable(contentWidth(), null);
    }

    /**
     * Add a paragraph of text at the current cursor position.
     */
    public void addParagraph(final String text, final ReportResources.TextStyle style) throws IOException
    {
        addParagraph(text, style.font(), style.fontSize(), style.color());
    }

    /**
     * Add a paragraph of text at the current cursor position.
     */
    public void addParagraph(final String text, final PDFont font, float fontSize, final Color color) throws IOException
    {
        ensurePage();

        float lineHeight = fontSize * 1.4f;

        if(!hasSpaceFor(lineHeight))
        {
            newPage();
        }

        try(PDPageContentStream cs = new PDPageContentStream(mDocument, mCurrentPage, PDPageContentStream.AppendMode.APPEND, true, true))
        {
            cs.beginText();
            cs.setFont(font, fontSize);
            cs.setNonStrokingColor(color);
            cs.newLineAtOffset(mMarginLeft, mCursorY - fontSize);
            cs.showText(text);
            cs.endText();
        }

        mCursorY -= lineHeight;
    }

    /**
     * Add the QC fail notice text (equivalent to old addQcFailNotice).
     */
    public void addQcFailNotice(final ReportResources reportResources) throws IOException
    {
        addParagraph(ReportResources.NOT_AVAILABLE, reportResources.tableContentStyle());
    }

    /**
     * Add an image at the current cursor position, centered horizontally.
     * The image is scaled to fit within maxWidth and maxHeight.
     */
    public void addImage(final String imagePath, float maxWidth, float maxHeight) throws IOException
    {
        addImage(imagePath, maxWidth, maxHeight, true);
    }

    /**
     * Add an image at the current cursor position.
     */
    public void addImage(final String imagePath, float maxWidth, float maxHeight, boolean centerHorizontally) throws IOException
    {
        ensurePage();

        PDImageXObject image = PDImageXObject.createFromFile(imagePath, mDocument);

        float imgWidth = image.getWidth();
        float imgHeight = image.getHeight();

        // Scale to fit within maxWidth and maxHeight while maintaining aspect ratio
        float scale = Math.min(maxWidth / imgWidth, maxHeight / imgHeight);
        if(scale > 1.0f)
        {
            scale = 1.0f; // don't upscale
        }

        float drawWidth = imgWidth * scale;
        float drawHeight = imgHeight * scale;

        if(!hasSpaceFor(drawHeight))
        {
            newPage();
        }

        float xPos = centerHorizontally
                ? mMarginLeft + (contentWidth() - drawWidth) / 2
                : mMarginLeft;

        float yPos = mCursorY - drawHeight;

        try(PDPageContentStream cs = new PDPageContentStream(mDocument, mCurrentPage, PDPageContentStream.AppendMode.APPEND, true, true))
        {
            cs.drawImage(image, xPos, yPos, drawWidth, drawHeight);
        }

        mCursorY = yPos - 5; // small spacing after image
    }

    /**
     * Create a Boxable BaseTable at a specific X position (for side-by-side layout).
     */
    public BaseTable createTableAtX(float width, float xStart) throws IOException
    {
        ensurePage();
        float yStart = mCursorY;
        float yStartNewPage = contentStartY();
        float bottomMargin = contentEndY();

        BaseTable table = new BaseTable(yStart, yStartNewPage, bottomMargin, width, xStart, mDocument, mCurrentPage, true, true);
        return table;
    }

    /**
     * Draw a table without advancing the cursor (for side-by-side layout).
     * Returns the final Y of this table so the caller can track both.
     */
    public float addTableNoAdvance(final BaseTable table) throws IOException
    {
        ensurePage();
        float finalY = table.draw();
        mCurrentPage = table.getCurrentPage();
        return finalY;
    }

    /**
     * Manually set the cursor Y position.
     */
    public void setCursorY(float y)
    {
        mCursorY = y;
    }

    /**
     * Draw an image at a specific X,Y position without advancing the cursor.
     * Returns the actual drawn height (for the caller to track row heights).
     */
    public float addImageAt(final String imagePath, float x, float y, float maxWidth, float maxHeight) throws IOException
    {
        ensurePage();

        PDImageXObject image = PDImageXObject.createFromFile(imagePath, mDocument);

        float imgWidth = image.getWidth();
        float imgHeight = image.getHeight();

        float scale = Math.min(maxWidth / imgWidth, maxHeight / imgHeight);
        if(scale > 1.0f)
        {
            scale = 1.0f;
        }

        float drawWidth = imgWidth * scale;
        float drawHeight = imgHeight * scale;

        // Center the image within the maxWidth cell
        float xPos = x + (maxWidth - drawWidth) / 2;
        float yPos = y - drawHeight;

        try(PDPageContentStream cs = new PDPageContentStream(mDocument, mCurrentPage, PDPageContentStream.AppendMode.APPEND, true, true))
        {
            cs.drawImage(image, xPos, yPos, drawWidth, drawHeight);
        }

        return drawHeight;
    }

    /**
     * Add vertical spacing.
     */
    public void addSpacing(float height) throws IOException
    {
        ensurePage();

        if(!hasSpaceFor(height))
        {
            newPage();
        }
        else
        {
            mCursorY -= height;
        }
    }

    /**
     * Write footers on all pages (second pass).
     */
    public void writeFooters() throws IOException
    {
        if(mPageEventHandler != null)
        {
            mPageEventHandler.writeFooters(mDocument);
        }
    }
}
