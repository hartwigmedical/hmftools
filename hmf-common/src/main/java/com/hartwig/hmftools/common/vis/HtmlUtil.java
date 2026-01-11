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
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.Pair;
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
        String readName = firstRead.getReadName();
        ChrBaseRegion alignment = new ChrBaseRegion(
                firstRead.getReferenceName(), firstRead.getAlignmentStart(), firstRead.getAlignmentEnd());
        ChrBaseRegion mateAlignment = null;
        if(firstRead.getReadPairedFlag() && !firstRead.getMateUnmappedFlag())
        {
            String mateChromosome = firstRead.getMateReferenceName();
            int mateAlignmentStart = firstRead.getMateAlignmentStart();
            int mateAlignmentEnd = getMateAlignmentEnd(firstRead);
            mateAlignment = new ChrBaseRegion(mateChromosome, mateAlignmentStart, mateAlignmentEnd);
        }

        String cigarStr = firstRead.getCigarString();
        String mateCigarStr = firstRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(mateCigarStr == null)
            mateCigarStr = "missing";

        int insertSize = inferredInsertSizeAbs(firstRead);
        String orientationStr = getOrientationString(firstRead);
        int mapQ = firstRead.getMappingQuality();
        int readNM = getReadNM(firstRead);
        String consensusTypeAttribute = firstRead.getStringAttribute(CONSENSUS_TYPE_ATTRIBUTE);
        String consensusReadAttribute = firstRead.hasAttribute(CONSENSUS_READ_ATTRIBUTE)
	        ? firstRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE)
	        : null;
        Integer secondMapQ = null;
        Integer secondReadNM = null;
        if(secondRead != null)
        {
            secondMapQ = secondRead.getMappingQuality();
            secondReadNM = getReadNM(secondRead);
        }

        return renderReadInfoTable(readName, alignment, mateAlignment, cigarStr, mateCigarStr, insertSize, orientationStr, mapQ, readNM, consensusTypeAttribute, consensusReadAttribute, secondMapQ, secondReadNM, null);
    }

    public static DomContent renderReadInfoTable(final String readName, final ChrBaseRegion alignment,
            @Nullable final ChrBaseRegion mateAlignment, final String cigarStr, @Nullable final String mateCigarStr,
            @Nullable final Integer insertSize, @Nullable final String orientationStr, @Nullable final Integer mapQ,
            @Nullable final Integer readNM, @Nullable final String consensusTypeAttribute, @Nullable final String consensusReadAttribute,
            @Nullable final Integer secondMapQ, @Nullable final Integer secondReadNM, @Nullable final List<Pair<String, String>> extraInfo)
    {
        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);
        CssBuilder readInfoStyle = baseDivStyle.display("none");

        List<DomContent> readInfoRows = Lists.newArrayList();

        readInfoRows.add(tr(td("Read name:"), td(readName)));
        String alignmentStr = format("%s:%s-%s", alignment.Chromosome, alignment.start(), alignment.end());

        String mateAlignmentStr = "unmapped";
        if(mateAlignment != null)
        {
            String mateChromosome = mateAlignment.Chromosome;
            int mateAlignmentStart = mateAlignment.start();
            int mateAlignmentEnd = mateAlignment.end();
            String mateAlignmentEndStr = mateAlignmentEnd == NO_POSITION ? "?" : String.valueOf(mateAlignmentEnd);
            mateAlignmentStr = format("%s:%d-%s", mateChromosome, mateAlignmentStart, mateAlignmentEndStr);
        }

        readInfoRows.add(tr(td("Alignment:"), td(alignmentStr + ", " + mateAlignmentStr)));

        readInfoRows.add(tr(td("Cigar:"), td(cigarStr + ", " + (mateCigarStr == null ? "missing" : mateCigarStr))));
        readInfoRows.add(tr(td("Insert size:"), td(insertSize == null ? "unknown" : String.valueOf(insertSize))));
        readInfoRows.add(tr(td("Orientation:"), td(orientationStr == null ? "unknown" : orientationStr)));

        String firstMapQStr = mapQ == null ? "unknown" : String.valueOf(mapQ);
        String secondMapQStr = secondMapQ == null ? "" : java.lang.String.valueOf(secondMapQ);
        String mapQStr = secondMapQ == null ? firstMapQStr : firstMapQStr + ", " + secondMapQStr;
        readInfoRows.add(tr(td("MapQ:"), td(mapQStr)));

        String firstNumMutationsStr = readNM == null ? "unknown" : String.valueOf(readNM);
        String secondNumMutationsStr = secondReadNM == null ? "" : java.lang.String.valueOf(secondReadNM);
        String numMutationsStr = secondReadNM == null ? firstNumMutationsStr : firstNumMutationsStr + ", " + secondNumMutationsStr;
        readInfoRows.add(tr(td("NM:"), td(numMutationsStr)));

        String umiTypeStr = consensusTypeAttribute;
        if(umiTypeStr != null)
            readInfoRows.add(tr(td("Dup type:"), td(umiTypeStr)));

        String dupCountStr = "0";
        if(consensusReadAttribute != null)
            dupCountStr = consensusReadAttribute.split(CONSENSUS_INFO_DELIM, 2)[0];

        readInfoRows.add(tr(td("Dup count:"), td(dupCountStr)));

        if(extraInfo != null)
        {
            for(Pair<String, String> entry : extraInfo)
                readInfoRows.add(tr(td(entry.getLeft()), td(entry.getRight())));
        }

        DomContent readInfoTable = styledTable(readInfoRows, CssBuilder.EMPTY);
        return div(readInfoTable).withClass("read-info").withStyle(readInfoStyle.toString());
    }
}
