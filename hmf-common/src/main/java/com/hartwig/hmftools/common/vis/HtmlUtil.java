package com.hartwig.hmftools.common.vis;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_INFO_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getNumEvents;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getOrientationString;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;

import static j2html.TagCreator.div;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.script;
import static j2html.TagCreator.table;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import j2html.tags.DomContent;

public final class HtmlUtil
{
    private HtmlUtil() {}

    // sizes
    private static final int BASE_FONT_SIZE = 10;

    // styles
    public static final CssBuilder BASE_FONT_STYLE = CssBuilder.EMPTY.fontSizePt(BASE_FONT_SIZE).fontFamily("sans-serif");

    public static final DomContent JQUERY_SCRIPT =
            rawHtml("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js\"></script>");

    private static final AtomicReference<DomContent> JAVASCRIPT = new AtomicReference<>(null);

    public static DomContent styledTable(final List<DomContent> elems, final CssBuilder style)
    {
        return table().with(elems).withStyle(BASE_FONT_STYLE.merge(style).toString());
    }

    public static DomContent getJavascript()
    {
        return JAVASCRIPT.updateAndGet((final DomContent currentRef) ->
        {
            if(currentRef != null)
                return currentRef;

            InputStream inputStream = HtmlUtil.class.getResourceAsStream("/vis/sagevis.js");
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            String scriptContent = reader.lines().collect(Collectors.joining("\n"));
            return script(rawHtml(scriptContent)).attr("type", "text/javascript");
        });
    }

    private static int getReadNM(final SAMRecord read) { return getNumEvents(read); }

    public static DomContent renderReadInfoTable(final SAMRecord firstRead, @Nullable final SAMRecord secondRead)
    {
        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);
        CssBuilder readInfoStyle = baseDivStyle.display("none");

        List<DomContent> readInfoRows = Lists.newArrayList();

        readInfoRows.add(tr(td("Read name:"), td(firstRead.getReadName())));

        String alignmentStr = format("%s:%s-%s", firstRead.getReferenceName(), firstRead.getAlignmentStart(), firstRead.getAlignmentEnd());
        String mateAlignmentStr = "unmapped";
        if(firstRead.getReadPairedFlag() && !firstRead.getMateUnmappedFlag())
        {
            String mateChromosome = firstRead.getMateReferenceName();
            int mateAlignmentStart = firstRead.getMateAlignmentStart();
            int mateAlignmentEnd = getMateAlignmentEnd(firstRead);
            String mateAlignmentEndStr = mateAlignmentEnd == NO_POSITION ? "?" : String.valueOf(mateAlignmentEnd);
            mateAlignmentStr = format("%s:%d-%s", mateChromosome, mateAlignmentStart, mateAlignmentEndStr);
        }
        readInfoRows.add(tr(td("Alignment:"), td(alignmentStr + ", " + mateAlignmentStr)));

        String mateCigarStr = firstRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(mateCigarStr == null)
        {
            mateCigarStr = "missing";
        }

        readInfoRows.add(tr(td("Cigar:"), td(firstRead.getCigarString() + ", " + mateCigarStr)));
        int insertSize = inferredInsertSizeAbs(firstRead);
        readInfoRows.add(tr(td("Insert size:"), td(String.valueOf(insertSize))));
        readInfoRows.add(tr(td("Orientation:"), td(getOrientationString(firstRead))));

        String firstMapQStr = String.valueOf(firstRead.getMappingQuality());
        String secondMapQStr = secondRead == null ? "" : String.valueOf(secondRead.getMappingQuality());
        String mapQStr = secondRead == null ? firstMapQStr : firstMapQStr + ", " + secondMapQStr;
        readInfoRows.add(tr(td("MapQ:"), td(mapQStr)));

        String firstNumMutationsStr = String.valueOf(getReadNM(firstRead));
        String secondNumMutationsStr = secondRead == null ? "" : String.valueOf(getReadNM(secondRead));
        String numMutationsStr = secondRead == null ? firstNumMutationsStr : firstNumMutationsStr + ", " + secondNumMutationsStr;
        readInfoRows.add(tr(td("NM:"), td(numMutationsStr)));

        String umiTypeStr = firstRead.getStringAttribute(CONSENSUS_TYPE_ATTRIBUTE);
        if(umiTypeStr != null)
        {
            readInfoRows.add(tr(td("Dup type:"), td(umiTypeStr)));
        }

        String dupCountStr = "0";
        if(firstRead.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
        {
            dupCountStr = firstRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE).split(CONSENSUS_INFO_DELIM, 2)[0];
        }

        readInfoRows.add(tr(td("Dup count:"), td(dupCountStr)));

        DomContent readInfoTable = styledTable(readInfoRows, CssBuilder.EMPTY);
        return div(readInfoTable).withClass("read-info").withStyle(readInfoStyle.toString());
    }
}
