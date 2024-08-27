package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.esvee.utils.vcfcompare.CompareConfig;
import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatch;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.Nullable;

public class LineLinkWriter
{
    public final String mSampleId;
    public final String mOutputDir;
    public final String mOutputId;

    private final boolean mIncludeNonPass;

    public LineLinkWriter(String sampleId, String outputDir, String outputId, boolean includeNonPass)
    {
        mSampleId = sampleId;
        mOutputDir = outputDir;
        mOutputId = outputId;

        mIncludeNonPass = includeNonPass;
    }

    public LineLinkWriter(CompareConfig config)
    {
        mSampleId = config.SampleId;
        mOutputDir = config.OutputDir;
        mOutputId = config.OutputId;

        mIncludeNonPass = config.IncludeNonPass;
    }

    private static final String OLD_PREFIX = "Old";
    private static final String NEW_PREFIX = "New";

    private static final List<String> FIXED_HEADER_FIELDS = List.of(
            "SampleId",
            "UnifiedPolyACoords",
            "PolyAMatchType",
            "HasAnyPass"
    );

    private static final List<String> HEADER_SUFFIXES = List.of(
            "VcfType",
            "PolyAHasLink",
            "PolyAId",
            "PolyARemoteId",
            "OtherId",
            "OtherRemoteId",
            "PolyACoords",
            "PolyARemoteCoords",
            "OtherCoords",
            "OtherRemoteCoords",
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
            SV_LOGGER.error("Failed to initialise output file: {}", e.toString());
            System.exit(1);
            return null;
        }
    }

    private class OutputRow
    {
        String VcfType = "";
        String PolyAHasLink = "";

        String PolyAId = "";
        String OtherId = "";
        String PolyARemoteId = "";
        String OtherRemoteId = "";

        String PolyACoords = "";
        String OtherCoords = "";
        String PolyARemoteCoords = "";
        String OtherRemoteCoords = "";

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

        public OutputRow(@Nullable VariantBreakend breakend)
        {
            if(breakend == null)
            {
                return;
            }

            boolean hasLink = breakend.hasLineLink();
            PolyAHasLink = String.valueOf(hasLink);

            LineLink lineLink;
            VariantBreakend polyASite;
            VariantBreakend otherSite;

            if(hasLink)
            {
                lineLink = breakend.LinkedLineBreakends;
                polyASite = lineLink.mPolyASite;
                otherSite = lineLink.mOtherSite;
            }
            else
            {
                lineLink = null;
                polyASite = breakend;
                otherSite = null;
            }

            VcfType = polyASite.SourceVcfType.toString();

            PolyAId = polyASite.Id;
            PolyACoords = polyASite.coordStr();
            PolyAInsertSeq = polyASite.InsertSequence;
            PolyASvType = polyASite.SvType;
            PolyAFilter = polyASite.filtersStr();
            PolyAQual = polyASite.qualStr();
            PolyAFrags = polyASite.fragsStr(mSampleId);

            if(!hasLink)
            {
                return;
            }

            if(lineLink.polyAHasRemote())
            {
                PolyARemoteCoords = polyASite.otherCoordStr();
                PolyARemoteId = polyASite.mateId();
            }

            OtherId = otherSite.Id;
            OtherCoords = otherSite.coordStr();
            OtherInsertSeq = otherSite.InsertSequence;
            OtherSvType = otherSite.SvType;
            OtherFilter = otherSite.filtersStr();
            OtherQual = otherSite.qualStr();
            OtherFrags = otherSite.fragsStr(mSampleId);

            if(lineLink.otherHasRemote())
            {
                OtherRemoteCoords = otherSite.otherCoordStr();
                OtherRemoteId = otherSite.mateId();
            }
        }

        public List<String> getOutputStrings()
        {
            return List.of(
                    VcfType,
                    PolyAHasLink,
                    PolyAId,
                    PolyARemoteId,
                    OtherId,
                    OtherRemoteId,
                    PolyACoords,
                    PolyARemoteCoords,
                    OtherCoords,
                    OtherRemoteCoords,
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
    }

    private boolean isLineInsertSiteOfInterest(@Nullable VariantBreakend breakend)
    {
        return breakend != null &&
                (mIncludeNonPass || breakend.isPassVariant()) &&
                breakend.isLineInsertionSite();
    }

    private static String getUnifiedPolyACoords(VariantBreakend oldBreakend, VariantBreakend newBreakend)
    {
        if(oldBreakend != null)
        {
            if(!oldBreakend.hasLineLink())
                return oldBreakend.coordStr();
            else
                return oldBreakend.LinkedLineBreakends.mPolyASite.coordStr();
        }
        else
        {
            if(!newBreakend.hasLineLink())
                return newBreakend.coordStr();
            else
                return newBreakend.LinkedLineBreakends.mPolyASite.coordStr();
        }
    }

    private static boolean variantHasAnyPass(@Nullable VariantBreakend breakend)
    {
        boolean anyPass;

        if(breakend == null)
        {
            anyPass = false;
        }
        else if(!breakend.hasLineLink())
        {
            anyPass = breakend.isPassVariant();
        }
        else
        {
            anyPass = breakend.LinkedLineBreakends.mPolyASite.isPassVariant();
            anyPass |= breakend.LinkedLineBreakends.mOtherSite.isPassVariant();
        }

        return anyPass;
    }

    public void writeBreakends(BreakendMatcher breakendMatcher)
    {
        try
        {
            BufferedWriter writer = initialiseWriter();

            List<BreakendMatch> breakendMatches = breakendMatcher.getBreakendMatches();

            for(BreakendMatch match : breakendMatches)
            {
                VariantBreakend oldBreakend = match.OldBreakend;
                VariantBreakend newBreakend = match.NewBreakend;

                if(!isLineInsertSiteOfInterest(oldBreakend) && !isLineInsertSiteOfInterest(newBreakend))
                    continue;

                List<String> rowStrings = new ArrayList<>();

                rowStrings.add(mSampleId);

                String unifiedPolyACoords = getUnifiedPolyACoords(oldBreakend, newBreakend);
                rowStrings.add(unifiedPolyACoords);

                String polyAMatchType = match.Type.toString();
                rowStrings.add(polyAMatchType);

                String hasAnyPass = String.valueOf(variantHasAnyPass(oldBreakend) || variantHasAnyPass(newBreakend));
                rowStrings.add(hasAnyPass);

                List<String> oldStrings = new OutputRow(oldBreakend).getOutputStrings();
                List<String> newStrings = new OutputRow(newBreakend).getOutputStrings();

                for(int i = 0; i < oldStrings.size(); i++)
                {
                    rowStrings.add(oldStrings.get(i));
                    rowStrings.add(newStrings.get(i));
                }

                writer.write(String.join(TSV_DELIM, rowStrings));
                writer.newLine();
            }

            FileWriterUtils.closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("Failed to write output file: {}", e.toString());
        }
    }
}
