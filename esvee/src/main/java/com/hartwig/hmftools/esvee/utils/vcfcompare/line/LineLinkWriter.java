package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatch;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.Nullable;

public class LineLinkWriter
{
    public final String mSampleId;
    public final String mOutputDir;
    public final String mOutputId;

    private final boolean mShowNonPass;

    private final BufferedWriter mWriter;

    public LineLinkWriter(String sampleId, String outputDir, String outputId, boolean showNonPass)
    {
        mSampleId = sampleId;
        mOutputDir = outputDir;
        mOutputId = outputId;

        mShowNonPass = showNonPass;

        mWriter = initialiseWriter();
    }

    private static final String OLD_PREFIX = "Old";
    private static final String NEW_PREFIX = "New";

    private static final List<String> FIXED_HEADER_FIELDS = List.of(
            "UnifiedPolyACoords",
            "PolyAMatchType"
    );

    private static final List<String> HEADER_SUFFIXES = List.of(
            "LineCoords",
            "PolyAId",
            "OtherId",
            "PolyACoords",
            "OtherCoords",
            "RemoteCoords",
            "PolyAInsertSeq",
            "OtherInsertSeq",
            "PolyASvType",
            "OtherSvType",
            "PolyAFilter",
            "OtherFilter",
            "PolyAQual",
            "OtherQual",
            "PolyAFrags",
            "OtherFrags"
    );

    private static String header()
    {
        List<String> headerStrings = new ArrayList<>();

        headerStrings.addAll(FIXED_HEADER_FIELDS);

        // Alternating
        for(String suffix : HEADER_SUFFIXES)
        {
            headerStrings.add(OLD_PREFIX + suffix);
            headerStrings.add(NEW_PREFIX + suffix);
        }

        return String.join(TSV_DELIM, headerStrings);
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".sv_compare.line";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += TSV_EXTENSION;

            SV_LOGGER.info("Writing LINE comparison file: {}", fileName);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(fileName, false);
            writer.write(header());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    public void closeWriter() { FileWriterUtils.closeBufferedWriter(mWriter); }

    private List<String> getOutputStrings(@Nullable VariantBreakend insertSite)
    {
        String LineCoords = "";

        String PolyAId = "";
        String OtherId = "";
        String PolyACoords = "";
        String OtherCoords = "";
        String RemoteCoords = "";
        String PolyAInsertSeq = "";
        String OtherInsertSeq = "";
        String PolyASvType = "";
        String OtherSvType = "";
        String PolyAFilter = "";
        String OtherFilter = "";
        String PolyAQual = "";
        String OtherQual = "";
        String PolyAFrags = "";
        String OtherFrags = "";

        //System.out.println(insertSite + "    " + insertSite.LinkedLineBreakend);

        if(insertSite != null)
        {
            VariantBreakend polyASite = insertSite;
            VariantBreakend otherSite = insertSite.LinkedLineBreakend;

            LineCoords = polyASite.lineLinkCoordStr();

            PolyAId = polyASite.id();
            PolyACoords = polyASite.coordStr();
            PolyAInsertSeq = polyASite.InsertSequence;
            PolyASvType = polyASite.SvType;
            PolyAFilter = polyASite.filtersStr();
            PolyAQual = polyASite.qualStr();
            PolyAFrags = polyASite.fragsStr(mSampleId);

            if(otherSite != null)
            {
                OtherId = otherSite.id();
                OtherCoords = otherSite.coordStr();
                OtherInsertSeq = otherSite.InsertSequence;

                if(otherSite.isTranslocation())
                    RemoteCoords = otherSite.otherCoordStr();

                OtherSvType = otherSite.SvType;
                OtherFilter = otherSite.filtersStr();
                OtherQual = otherSite.qualStr();
                OtherFrags = otherSite.fragsStr(mSampleId);
            }
        }

        return List.of(
                LineCoords,
                PolyAId,
                OtherId,
                PolyACoords,
                OtherCoords,
                RemoteCoords,
                PolyAInsertSeq,
                OtherInsertSeq,
                PolyASvType,
                OtherSvType,
                PolyAFilter,
                OtherFilter,
                PolyAQual,
                OtherQual,
                PolyAFrags,
                OtherFrags
        );
    }

    private boolean isLineInsertSiteOfInterest(@Nullable VariantBreakend breakend)
    {
        return breakend != null &&
                (mShowNonPass || breakend.isPassVariant()) &&
                breakend.isLineInsertionSite();
    }

    public void writeBreakends(BreakendMatcher breakendMatcher)
    {
        try
        {
            SV_LOGGER.info("Writing LINE breakends");

            List<BreakendMatch> breakendMatches = breakendMatcher.getBreakendMatches();

            for(BreakendMatch match : breakendMatches)
            {
                VariantBreakend oldBreakend = match.OldBreakend;
                VariantBreakend newBreakend = match.NewBreakend;

                // System.out.println(oldBreakend + "    " + newBreakend);
                // int a = 0;

                if(!isLineInsertSiteOfInterest(oldBreakend) && !isLineInsertSiteOfInterest(newBreakend))
                    continue;

                List<String> rowStrings = new ArrayList<>();

                String unifiedPolyACoords = oldBreakend != null ? oldBreakend.coordStr() : newBreakend.coordStr();

                String polyAMatchType = match.Type.toString();

                rowStrings.add(unifiedPolyACoords);
                rowStrings.add(polyAMatchType);

                List<String> oldStrings = getOutputStrings(oldBreakend);
                List<String> newStrings = getOutputStrings(newBreakend);

                for(int i = 0; i < oldStrings.size(); i++)
                {
                    rowStrings.add(oldStrings.get(i));
                    rowStrings.add(newStrings.get(i));
                }

                mWriter.write(String.join(TSV_DELIM, rowStrings));
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("Failed to write output file: {}", e.toString());
        }
    }
}
