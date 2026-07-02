package com.hartwig.hmftools.orange.report;

import java.awt.Color;
import java.io.IOException;
import java.io.InputStream;

import com.hartwig.hmftools.orange.OrangeApplication;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.font.PDFont;
import org.apache.pdfbox.pdmodel.font.PDType0Font;

import be.quodlibet.boxable.Cell;

public class ReportResources
{
    static final String METADATA_TITLE = "HMF ORANGE Report v" + OrangeApplication.VERSION;
    static final String METADATA_AUTHOR = "Hartwig Pipeline";

    public static final String NOT_AVAILABLE = "NA";
    public static final String NONE = "NONE";

    public static final float PAGE_MARGIN_TOP = 100; // Top margin also excludes the chapter title, which is rendered in the header
    public static final float PAGE_MARGIN_LEFT = 30;
    public static final float PAGE_MARGIN_RIGHT = 30;
    public static final float PAGE_MARGIN_BOTTOM = 40;

    public static final int HEADER_ORANGE_HEIGHT = 60;
    public static final int FOOTER_HEIGHT = 20;

    public static final int FRONT_CIRCOS_IMAGE_HEIGHT = 400;
    public static final int FULL_PAGE_IMAGE_WIDTH = 750;
    public static final int FULL_PAGE_IMAGE_HEIGHT = 430;

    public static final Color PALETTE_WHITE = new Color(255, 255, 255);
    public static final Color PALETTE_BLACK = new Color(0, 0, 0);

    public static final Color PALETTE_DARK_GREY = new Color(39, 47, 50);
    public static final Color PALETTE_MID_GREY = new Color(101, 106, 108);
    public static final Color PALETTE_LIGHT_GREY = new Color(211, 211, 211);
    public static final Color PALETTE_GAINSBORO_GREY = new Color(220, 220, 220);
    private static final int GREY_FACTOR = 240;
    public static final Color PALETTE_SMOKE_GREY = new Color(GREY_FACTOR, GREY_FACTOR, GREY_FACTOR);
    public static final Color PALETTE_BLUE = new Color(38, 90, 166);

    public static final Color PALETTE_ORANGE = new Color(242, 139, 31);
    public static final Color PALETTE_ORANGE_1 = new Color(255, 165, 0);
    public static final Color PALETTE_ORANGE_2 = new Color(235, 155, 0);
    public static final Color PALETTE_ORANGE_3 = new Color(215, 145, 0);
    public static final Color PALETTE_ORANGE_4 = new Color(195, 135, 0);
    public static final Color PALETTE_ORANGE_5 = new Color(175, 125, 0);
    public static final Color PALETTE_ORANGE_6 = new Color(155, 115, 0);

    private static final String FONT_REGULAR_PATH = "fonts/nimbus-sans/NimbusSansL-Regular.ttf";
    private static final String FONT_BOLD_PATH = "fonts/nimbus-sans/NimbusSansL-Bold.ttf";

    // Standard font sizes
    public static final float FONT_SIZE_CHAPTER_TITLE = 10;
    public static final float FONT_SIZE_TABLE_TITLE = 8;
    public static final float FONT_SIZE_TABLE_HEADER = 7;
    public static final float FONT_SIZE_TABLE_CONTENT = 7;
    public static final float FONT_SIZE_KEY = 7;
    public static final float FONT_SIZE_VALUE = 7;
    public static final float FONT_SIZE_SUB_TEXT = 6;
    public static final float FONT_SIZE_PAGE_NUMBER = 7;
    public static final float FONT_SIZE_DISCLAIMER = 6;
    public static final float FONT_SIZE_QC_WARNING = 7;
    public static final float FONT_SIZE_SIDE_PANEL_LABEL = 7;
    public static final float FONT_SIZE_SIDE_PANEL_VALUE = 10;
    public static final float FONT_SIZE_URL = 7;

    private final PDFont mFontRegular;
    private final PDFont mFontBold;

    private ReportResources(final PDFont fontRegular, final PDFont fontBold)
    {
        mFontRegular = fontRegular;
        mFontBold = fontBold;
    }

    public static ReportResources create(final PDDocument document)
    {
        return new ReportResources(loadFont(document, FONT_REGULAR_PATH), loadFont(document, FONT_BOLD_PATH));
    }

    public PDFont fontRegular()
    {
        return mFontRegular;
    }

    public PDFont fontBold()
    {
        return mFontBold;
    }

    // Style accessors returning TextStyle
    public TextStyle chapterTitleStyle() { return new TextStyle(mFontBold, FONT_SIZE_CHAPTER_TITLE, PALETTE_ORANGE); }
    public TextStyle tableTitleStyle() { return new TextStyle(mFontBold, FONT_SIZE_TABLE_TITLE, PALETTE_ORANGE); }
    public TextStyle tableHeaderStyle() { return new TextStyle(mFontRegular, FONT_SIZE_TABLE_HEADER, PALETTE_MID_GREY); }
    public TextStyle tableContentStyle() { return new TextStyle(mFontRegular, FONT_SIZE_TABLE_CONTENT, PALETTE_DARK_GREY); }
    public TextStyle keyStyle() { return new TextStyle(mFontRegular, FONT_SIZE_KEY, PALETTE_MID_GREY); }
    public TextStyle valueStyle() { return new TextStyle(mFontRegular, FONT_SIZE_VALUE, PALETTE_MID_GREY); }
    public TextStyle subTextStyle() { return new TextStyle(mFontRegular, FONT_SIZE_SUB_TEXT, PALETTE_BLACK); }
    public TextStyle pageNumberStyle() { return new TextStyle(mFontBold, FONT_SIZE_PAGE_NUMBER, PALETTE_ORANGE); }
    public TextStyle disclaimerStyle() { return new TextStyle(mFontRegular, FONT_SIZE_DISCLAIMER, PALETTE_MID_GREY); }
    public TextStyle qcWarningStyle() { return new TextStyle(mFontBold, FONT_SIZE_QC_WARNING, PALETTE_DARK_GREY); }
    public TextStyle sidePanelLabelStyle() { return new TextStyle(mFontBold, FONT_SIZE_SIDE_PANEL_LABEL, PALETTE_WHITE); }
    public TextStyle sidePanelValueStyle() { return new TextStyle(mFontBold, FONT_SIZE_SIDE_PANEL_VALUE, PALETTE_WHITE); }
    public TextStyle urlStyle() { return new TextStyle(mFontRegular, FONT_SIZE_URL, PALETTE_BLUE); }

    @SuppressWarnings("unchecked")
    public void shadeCandidateCells(final java.util.List<? extends Cell> cells)
    {
        cells.forEach(x -> x.setFillColor(PALETTE_SMOKE_GREY));
    }

    private static PDFont loadFont(final PDDocument document, final String resourcePath)
    {
        try(InputStream is = ReportResources.class.getClassLoader().getResourceAsStream(resourcePath))
        {
            if(is == null)
            {
                throw new IllegalStateException("Font resource not found: " + resourcePath);
            }
            return PDType0Font.load(document, is);
        }
        catch(IOException exception)
        {
            throw new IllegalStateException(exception);
        }
    }

    /** Simple holder for font + size + color styling. */
    public static class TextStyle
    {
        private final PDFont mFont;
        private final float mFontSize;
        private final Color mColor;

        public TextStyle(final PDFont font, float fontSize, final Color color)
        {
            mFont = font;
            mFontSize = fontSize;
            mColor = color;
        }

        public PDFont font() { return mFont; }
        public float fontSize() { return mFontSize; }
        public Color color() { return mColor; }
    }
}
