package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.esvee.utils.vcfcompare.CompareConfig;
import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatch;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LineLinkWriter
{
    public final String mSampleId;
    public final String mOutputDir;
    public final String mOutputId;

    private final boolean mIncludeNonPass;

    private final BreakendMatcher mBreakendMatcher;

    public LineLinkWriter(BreakendMatcher breakendMatcher, String sampleId, String outputDir, String outputId, boolean includeNonPass)
    {
        mSampleId = sampleId;
        mOutputDir = outputDir;
        mOutputId = outputId;

        mIncludeNonPass = includeNonPass;

        mBreakendMatcher = breakendMatcher;
    }

    public LineLinkWriter(BreakendMatcher breakendMatcher, CompareConfig config)
    {
        mSampleId = config.SampleId;
        mOutputDir = config.OutputDir;
        mOutputId = config.OutputId;

        mIncludeNonPass = config.IncludeNonPass;

        mBreakendMatcher = breakendMatcher;
    }

    private static final String OLD_PREFIX = "Old";
    private static final String NEW_PREFIX = "New";

    private static final List<String> FIXED_HEADER_FIELDS = List.of(
            "SampleId",
            "UnifiedPolyACoords",
            "UnifiedPolyACoordsCount",
            "PolyAMatchType",
            "HasAnyPass"
    );

    private static final List<String> HEADER_SUFFIXES = List.of(
            "VcfType",
            "LinkType",
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

    private class BreakendRowValues
    {
        String VcfType = "";
        String LinkType = "";

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

        private boolean mHasLink = false;
        private LineLink mLineLink = null;

        VariantBreakend mPolyASite;
        VariantBreakend mOtherSite;

        public BreakendRowValues(@Nullable VariantBreakend breakend)
        {
            if(breakend == null)
                return;

            mHasLink = breakend.hasLineLink();

            if(mHasLink)
            {
                LinkType = breakend.LinkedLineBreakends.mType.toString();

                mLineLink = breakend.LinkedLineBreakends;
                mPolyASite = mLineLink.mPolyASite;
                mOtherSite = mLineLink.mOtherSite;
            }
            else
            {
                LinkType = LineLinkType.NO_LINK.toString();

                mLineLink = null;
                mPolyASite = breakend;
                mOtherSite = null;
            }

            setPolyASiteValues(mPolyASite);

            if(mOtherSite != null)
                setOtherSiteValues(mOtherSite);
        }

        private void setPolyASiteValues(@NotNull VariantBreakend polyASite)
        {
            VcfType = polyASite.SourceVcfType.toString();

            PolyAId = polyASite.Id;
            PolyACoords = polyASite.coordStr();
            PolyAInsertSeq = polyASite.InsertSequence;
            PolyASvType = polyASite.SvType;
            PolyAFilter = polyASite.filtersStr();
            PolyAQual = polyASite.qualStr();
            PolyAFrags = polyASite.fragsStr(mSampleId);

            if(mHasLink && mLineLink.polyAHasRemote())
            {
                PolyARemoteCoords = polyASite.otherCoordStr();
                PolyARemoteId = polyASite.mateId();
            }
        }

        private void setOtherSiteValues(@NotNull VariantBreakend otherSite)
        {
            OtherId = otherSite.Id;
            OtherCoords = otherSite.coordStr();
            OtherInsertSeq = otherSite.InsertSequence;
            OtherSvType = otherSite.SvType;
            OtherFilter = otherSite.filtersStr();
            OtherQual = otherSite.qualStr();
            OtherFrags = otherSite.fragsStr(mSampleId);

            if(mHasLink && mLineLink.otherHasRemote())
            {
                OtherRemoteCoords = otherSite.otherCoordStr();
                OtherRemoteId = otherSite.mateId();
            }
        }

        private void setInferredOtherSiteValues(@NotNull VariantBreakend otherBreakend)
        {
            LineLink inferredLink = otherBreakend.InferredLinkedLineBreakends;

            SV_LOGGER.trace("Used inferred link in breakend[{}] to assign otherSite[{}]", otherBreakend, inferredLink.mOtherSite);
            setPolyASiteValues(otherBreakend);
            setOtherSiteValues(inferredLink.mOtherSite);

            LinkType = inferredLink.mType.toString();
        }

        public List<String> getOutputStrings()
        {
            return List.of(
                    VcfType,
                    LinkType,
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

    private class RowValues
    {
        BreakendRowValues mBreakend1RowValues;
        BreakendRowValues mBreakend2RowValues;

        public RowValues(@Nullable VariantBreakend breakend1, @Nullable VariantBreakend breakend2)
        {
            if(breakend1 == null && breakend2 == null)
                throw new IllegalStateException("`breakend1` and `breakend2` cannot both be null");

            mBreakend1RowValues = new BreakendRowValues(breakend1);
            mBreakend2RowValues = new BreakendRowValues(breakend2);

            if(breakend2 == null && breakend1.hasInferredLineLink())
            {
                mBreakend2RowValues.setInferredOtherSiteValues(breakend1);
            }

            if(breakend1 == null && breakend2.hasInferredLineLink())
            {
                mBreakend1RowValues.setInferredOtherSiteValues(breakend2);
            }
        }

        public List<String> getAlternatingStrings(){

            List<String> breakend1Strings = mBreakend1RowValues.getOutputStrings();
            List<String> breakend2Strings = mBreakend2RowValues.getOutputStrings();

            List<String> rowStrings = new ArrayList<>();
            for(int i = 0; i < breakend1Strings.size(); i++)
            {
                rowStrings.add(breakend1Strings.get(i));
                rowStrings.add(breakend2Strings.get(i));
            }

            return rowStrings;
        }
    }

    private boolean isLineInsertSiteOfInterest(@Nullable VariantBreakend breakend)
    {
        return breakend != null &&
                (mIncludeNonPass || breakend.isPassVariant()) &&
                breakend.hasPolyATail();
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

    private HashMap<String, Integer> countUnifiedPolyACoords(List<BreakendMatch> breakendMatches)
    {
        HashMap<String, Integer> countMap = new HashMap<>();

        for(BreakendMatch match : breakendMatches)
        {
            VariantBreakend oldBreakend = match.OldBreakend;
            VariantBreakend newBreakend = match.NewBreakend;

            if(!isLineInsertSiteOfInterest(oldBreakend) && !isLineInsertSiteOfInterest(newBreakend))
                continue;

            String coords = getUnifiedPolyACoords(oldBreakend, newBreakend);

            countMap.merge(coords, 1, Integer::sum);
        }

        return countMap;
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

    public void writeBreakends()
    {
        try
        {
            BufferedWriter writer = initialiseWriter();

            List<BreakendMatch> breakendMatches = mBreakendMatcher.getBreakendMatches();

            HashMap<String, Integer> unifiedPolyACoordsCountMap = countUnifiedPolyACoords(breakendMatches);

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

                Integer unifiedPolyACoordsCount = unifiedPolyACoordsCountMap.get(unifiedPolyACoords);
                rowStrings.add(unifiedPolyACoordsCount.toString());

                String polyAMatchType = match.Type.toString();
                rowStrings.add(polyAMatchType);

                String hasAnyPass = String.valueOf(variantHasAnyPass(oldBreakend) || variantHasAnyPass(newBreakend));
                rowStrings.add(hasAnyPass);

                RowValues rowValuesBreakendPair = new RowValues(oldBreakend, newBreakend);
                rowStrings.addAll(rowValuesBreakendPair.getAlternatingStrings());

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
