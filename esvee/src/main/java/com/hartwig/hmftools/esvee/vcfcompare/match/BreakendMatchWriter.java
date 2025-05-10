package com.hartwig.hmftools.esvee.vcfcompare.match;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.esvee.vcfcompare.CompareConfig;
import com.hartwig.hmftools.esvee.vcfcompare.VariantBreakend;
import com.hartwig.hmftools.esvee.vcfcompare.VcfType;

public class BreakendMatchWriter
{
    private final String mReferenceId;

    private final CompareConfig mConfig;

    private final List<BreakendMatch> mBreakendMatches;

    public BreakendMatchWriter(final List<BreakendMatch> breakendMatches, final CompareConfig config)
    {
        mBreakendMatches = breakendMatches;

        mConfig = config;
        mReferenceId = inferNormalId(config.SampleId);
    }

    private static String inferNormalId(String tumorId)
    {
        Pattern pattern = Pattern.compile("^[a-zA-Z0-9]+T");
        Matcher matcher = pattern.matcher(tumorId);

        String trimmedTumorId = null;
        if(matcher.find())
        {
            trimmedTumorId = matcher.group();
        }

        if(trimmedTumorId == null)
        {
            throw new IllegalStateException("Could not infer normal id from tumor id: " + tumorId);
        }

        return trimmedTumorId.substring(0, trimmedTumorId.length()-1) + "R";
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.formFilename("breakend");

            SV_LOGGER.debug("writing comparison file: {}", fileName);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(fileName, false);

            String header = String.join(
                    TSV_DELIM,
                    "OldId", "NewId", "MatchType", "Diffs", "OldSvCoords", "NewSvCoords", "OldCoords", "NewCoords",
                    "OldCipos", "NewCipos", "OldIhompos",    "NewIhompos", "OldHomSeq", "NewHomSeq",
                    "OldInsSeq", "NewInsSeq", "OldSvType", "NewSvType", "OldFilter", "NewFilter", "OldVcfType", "NewVcfType",
                    "OldQual", "NewQual", "OldTumorFrags", "NewTumorFrags", "OldNormalFrags","NewNormalFrags", "OldIsLine", "NewIsLine"
            );

            writer.write(header);
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

    public void writeBreakends()
    {
        BufferedWriter writer = initialiseWriter();

        for(BreakendMatch breakendMatch : mBreakendMatches)
        {
            if(!mConfig.WriteMatches && breakendMatch.Type == MatchType.EXACT_MATCH)
                continue;

            writeBreakend(writer, breakendMatch.OldBreakend, breakendMatch.NewBreakend, breakendMatch.Type);
        }

        FileWriterUtils.closeBufferedWriter(writer);
    }

    private void writeBreakend(
            final BufferedWriter writer, final VariantBreakend oldBreakend, final VariantBreakend newBreakend, final MatchType matchType)
    {
        String oldId = "";
        String oldSvCoords = "";
        String oldCoords = "";
        String oldCipos = "";
        String oldIhompos = "";
        String oldHomSeq = "";
        String oldInsSeq = "";
        String oldSvtype = "";
        String oldFilter = "";
        String oldVcfType = "";
        String oldQual = "";
        String oldTumorFrags  = "";
        String oldNormalFrags = "";
        String oldIsLine = "";
        if(oldBreakend != null)
        {
            oldId = oldBreakend.Context.getID();
            oldSvCoords = oldBreakend.svCoordStr();
            oldCoords = oldBreakend.coordStr();
            oldCipos = Arrays.toString(oldBreakend.Cipos);
            oldIhompos = Arrays.toString(oldBreakend.Ihompos);
            oldHomSeq = oldBreakend.Homseq;
            oldInsSeq = oldBreakend.InsertSequence;
            oldSvtype = oldBreakend.SvType;
            oldFilter = oldBreakend.filtersStr();
            oldVcfType = oldBreakend.SourceVcfType.toString();
            oldQual = oldBreakend.qualStr();
            oldTumorFrags = oldBreakend.fragsStr(mConfig.SampleId);
            oldNormalFrags = oldBreakend.fragsStr(mConfig.ReferenceId);
            oldIsLine = String.valueOf(oldBreakend.hasLineInfoFlag());
        }

        String newId = "";
        String newSvCoords = "";
        String newCoords = "";
        String newCipos = "";
        String newIhompos = "";
        String newHomSeq = "";
        String newInsSeq = "";
        String newSvtype = "";
        String newFilter = "";
        String newVcfType = "";
        String newQual = "";
        String newTumorFrags = "";
        String newNormalFrags = "";
        String newIsLine = "";

        if(newBreakend != null)
        {
            newId = newBreakend.Context.getID();
            newSvCoords = newBreakend.svCoordStr();
            newCoords = newBreakend.coordStr();
            newCipos = Arrays.toString(newBreakend.Cipos);
            newIhompos = Arrays.toString(newBreakend.Ihompos);
            newHomSeq = newBreakend.Homseq;
            newInsSeq = newBreakend.InsertSequence;
            newSvtype = newBreakend.SvType;
            newFilter = newBreakend.filtersStr();
            newVcfType = newBreakend.SourceVcfType.toString();
            newQual = newBreakend.qualStr();
            newTumorFrags = newBreakend.fragsStr(mConfig.SampleId);
            newNormalFrags = newBreakend.fragsStr(mReferenceId);
            newIsLine = String.valueOf(newBreakend.hasLineInfoFlag());
        }

        String diffs = "";
        if(oldBreakend != null && newBreakend != null)
        {
            diffs = String.join(",", getBreakendAttributeDiffs(oldBreakend, newBreakend, matchType));
        }

        try
        {
            String line = String.join(
                    TSV_DELIM,
                    oldId, newId,
                    matchType.toString(),
                    diffs, oldSvCoords, newSvCoords, oldCoords, newCoords,
                    oldCipos, newCipos, oldIhompos, newIhompos, oldHomSeq, newHomSeq, oldInsSeq, newInsSeq,
                    oldSvtype, newSvtype, oldFilter, newFilter, oldVcfType, newVcfType, oldQual, newQual,
                    oldTumorFrags, newTumorFrags, oldNormalFrags, newNormalFrags, oldIsLine, newIsLine
            );

            writer.write(line);
            writer.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("Failed to write output file: {}", e.toString());
        }
    }

    private static final String DIFF_PASS = "PASS_FILTER";
    private static final String DIFF_COORDS = "COORDS";
    private static final String DIFF_INSSEQ = "INSSEQ";
    private static final String DIFF_VCF_TYPE = "VCF_TYPE";
    private static final String DIFF_TUMOR_FRAGS = "VF_TUMOR";
    private static final String DIFF_NORMAL_FRAGS = "VF_NORMAL";

    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;

    private static boolean hasDiffWithinTolerance(double value1, double value2)
    {
        double diff = abs(value1 - value2);
        double diffPerc = diff / max(value1, value2);
        return diff > DEFAULT_MAX_DIFF && diffPerc > DEFAULT_MAX_DIFF_PERC;
    }

    // TODO: Make this a method in `BreakendMatch`
    private List<String> getBreakendAttributeDiffs(final VariantBreakend oldBreakend, final VariantBreakend newBreakend, final MatchType matchType)
    {
        List<String> diffSet = new ArrayList<>();

        if((oldBreakend.isPassVariant() & !newBreakend.isPassVariant()) || (!oldBreakend.isPassVariant() & newBreakend.isPassVariant()))
        {
            diffSet.add(DIFF_PASS);
        }

        // Coords
        if((matchType == MatchType.APPROX_MATCH || matchType == MatchType.COORDS_ONLY))
        {
            if(!oldBreakend.coordStr().equals(newBreakend.coordStr()))
                diffSet.add(DIFF_COORDS);

            if(!Arrays.equals(oldBreakend.Cipos, newBreakend.Cipos))
                diffSet.add(CIPOS);

            if(!Arrays.equals(oldBreakend.Ihompos, newBreakend.Ihompos))
                diffSet.add(IHOMPOS);

            if(!oldBreakend.Homseq.equals(newBreakend.Homseq))
                diffSet.add(HOMSEQ);

            if(!oldBreakend.InsertSequence.equals(newBreakend.InsertSequence))
                diffSet.add(DIFF_INSSEQ);

            if(!oldBreakend.SvType.equals(newBreakend.SvType))
                diffSet.add(SV_TYPE);
        }

        // VCF type
        VcfType oldSourceVcfType = (oldBreakend.SourceVcfType == VcfType.TRUTH) ? VcfType.SOMATIC : oldBreakend.SourceVcfType;
        VcfType newSourceVcfType = (newBreakend.SourceVcfType == VcfType.TRUTH) ? VcfType.SOMATIC : newBreakend.SourceVcfType;

        if(!oldSourceVcfType.equals(newSourceVcfType))
            diffSet.add(DIFF_VCF_TYPE);

        // VF
        double oldTumorFrags = oldBreakend.getExtendedAttributeAsDouble(mConfig.SampleId, TOTAL_FRAGS);
        double newTumorFrags = newBreakend.getExtendedAttributeAsDouble(mConfig.SampleId, TOTAL_FRAGS);

        if(hasDiffWithinTolerance(oldTumorFrags, newTumorFrags))
            diffSet.add(DIFF_TUMOR_FRAGS);

        double oldNormalFrags = oldBreakend.getExtendedAttributeAsDouble(mReferenceId, TOTAL_FRAGS);
        double newNormalFrags = newBreakend.getExtendedAttributeAsDouble(mReferenceId, TOTAL_FRAGS);

        if(hasDiffWithinTolerance(oldNormalFrags, newNormalFrags))
            diffSet.add(DIFF_NORMAL_FRAGS);

        return diffSet;
    }
}
