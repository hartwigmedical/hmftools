package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.sage.vis.ColorUtil.DARK_GREEN;
import static com.hartwig.hmftools.sage.vis.ColorUtil.PURPLE;
import static com.hartwig.hmftools.sage.vis.ColorUtil.interpolateColors;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.MAX_MAPQ_SHADING_CUTTOFF;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.READ_HEIGHT_PX;
import static com.hartwig.hmftools.sage.vis.SvgRender.FORWARD_STRAND_COLOR;
import static com.hartwig.hmftools.sage.vis.SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR;
import static com.hartwig.hmftools.sage.vis.SvgRender.REVERSE_STRAND_COLOR;
import static com.hartwig.hmftools.sage.vis.SvgRender.renderColoredBox;

import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.td;

import java.awt.Color;
import java.util.function.Function;

import htsjdk.samtools.SAMRecord;
import j2html.tags.specialized.TdTag;

public class ReadTableColumn
{
    public static final ReadTableColumn MATE_TYPE_COL = new ReadTableColumn("Mate", (final ReadEvidenceRecord record) ->
    {
        SAMRecord read = record.Fragment == null ? record.Read : record.Fragment.First;
        if(read.getMateUnmappedFlag())
        {
            return new ContentAndStyle(td("UN"), CssBuilder.EMPTY);
        }

        if(read.getReferenceName() != read.getMateReferenceName())
        {
            return new ContentAndStyle(td("TRL"), CssBuilder.EMPTY.color(PURPLE));
        }

        if(read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
        {
            return new ContentAndStyle(td("INV"), CssBuilder.EMPTY.color(Color.BLUE));
        }

        // dup
        boolean isDup = false;
        if(read.getReadNegativeStrandFlag())
        {
            isDup = read.getAlignmentEnd() < read.getMateAlignmentStart();
        }
        else
        {
            int mateAlignmentEnd = getMateAlignmentEnd(read);
            isDup = mateAlignmentEnd != NO_POSITION && read.getAlignmentStart() > mateAlignmentEnd;
        }

        if(isDup)
        {
            return new ContentAndStyle(td("DUP"), CssBuilder.EMPTY.color(DARK_GREEN));
        }

        if(abs(read.getInferredInsertSize()) >= SageVisConstants.INSERT_SIZE_CUTTOFF)
        {
            return new ContentAndStyle(td("DEL"), CssBuilder.EMPTY.color(Color.RED));
        }

        return new ContentAndStyle(td(), CssBuilder.EMPTY);
    });
    public static final ReadTableColumn MAP_QUAL_COL = new ReadTableColumn("MAPQ", (final ReadEvidenceRecord record) ->
    {
        int mapq = record.Read.getMappingQuality();
        double s = 1.0 * mapq / MAX_MAPQ_SHADING_CUTTOFF;
        Color color = interpolateColors(Color.RED, Color.GREEN, s);
        return new ContentAndStyle(td(rawHtml(renderColoredBox(READ_HEIGHT_PX, color).getSVGElement())), CssBuilder.EMPTY);
    });
    public static final ReadTableColumn FINAL_QUAL_COL = intValueColumn("FinalQ", (final ReadEvidenceRecord record) ->
    {
        return record.Qualities == null ? null : (int) Math.round(record.Qualities.ModifiedQuality);
    });
    public static final ReadTableColumn MOD_BASE_QUAL_COL = intValueColumn("modBQ", (final ReadEvidenceRecord record) ->
    {
        return record.Qualities == null ? null : (int) Math.round(record.Qualities.ModifiedBaseQuality);
    });
    public static final ReadTableColumn MOD_MAP_QUAL_COL = intValueColumn("modMQ", (final ReadEvidenceRecord record) ->
    {
        return record.Qualities == null ? null : record.Qualities.ModifiedMapQuality;
    });
    public static final ReadTableColumn RAW_BASE_QUAL_COL = intValueColumn("rawBQ", (final ReadEvidenceRecord record) ->
    {
        return record.Qualities == null ? null : (int) Math.round(record.Qualities.CalcBaseQuality);
    });
    public static final ReadTableColumn ORIENTATION_COL = new ReadTableColumn("Orientation", (final ReadEvidenceRecord record) ->
    {
        ReadEvidenceRecord.Orientation orientation = record.orientation();
        Color color = OVERLAPPING_FRAGMENT_BORDER_COLOR;
        switch(orientation)
        {
            case FORWARD:
                color = FORWARD_STRAND_COLOR;
                break;
            case REVERSE:
                color = REVERSE_STRAND_COLOR;
        }

        return new ContentAndStyle(td(rawHtml(renderColoredBox(READ_HEIGHT_PX, color).getSVGElement())), CssBuilder.EMPTY);
    });

    public final String Header;
    private final Function<ReadEvidenceRecord, ContentAndStyle> mGetter;

    private ReadTableColumn(final String header, final Function<ReadEvidenceRecord, ContentAndStyle> getter)
    {
        Header = header;
        mGetter = getter;
    }

    private static ReadTableColumn intValueColumn(final String header, final Function<ReadEvidenceRecord, Integer> valueGetter)
    {
        return new ReadTableColumn(header, (final ReadEvidenceRecord record) ->
        {
            CssBuilder style = CssBuilder.EMPTY.textAlign("right").paddingRight(CssSize.em(0.5f));
            Integer value = valueGetter.apply(record);
            if(value == null)
            {
                return new ContentAndStyle(td(), style);
            }
            else
            {
                return new ContentAndStyle(td(String.valueOf(value)), style);
            }
        });
    }

    public ContentAndStyle getContentAndStyle(final ReadEvidenceRecord readEvidence)
    {
        return mGetter.apply(readEvidence);
    }

    public static class ContentAndStyle
    {
        public final TdTag Content;
        public final CssBuilder Style;

        public ContentAndStyle(final TdTag content, final CssBuilder style)
        {
            Content = content;
            Style = style;
        }
    }
}
