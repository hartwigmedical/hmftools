package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.common.FileCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.vcfcompare.CoordMatchType.EXACT;
import static com.hartwig.hmftools.esvee.vcfcompare.CoordMatchType.NONE;
import static com.hartwig.hmftools.esvee.vcfcompare.CoordMatchType.determineMatchType;
import static com.hartwig.hmftools.esvee.vcfcompare.DiffUtils.checkValueDifference;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VcfComparer
{
    private final CompareConfig mConfig;

    private BufferedWriter mWriter;

    private final Set<Variant> mComparedVariants;

    public VcfComparer(final CompareConfig config)
    {
        mConfig = config;
        mWriter = null;
        mComparedVariants = Sets.newHashSet();
    }

    public void run()
    {
        // establish sample info from VCFs

        // SV_LOGGER.info("comparing VCFs for sample: " + mSampleId);
        VariantCache oldVariants = new VariantCache(mConfig.OldVcf);
        VariantCache newVariants = new VariantCache(mConfig.NewVcf);

        mWriter = initialiseWriter(oldVariants.genotypeIds().TumorId);

        VariantCache oldFilteredVariants = new VariantCache(mConfig.OldUnfilteredVcf);
        VariantCache newFilteredVariants = new VariantCache(mConfig.NewUnfilteredVcf);

        // begin matching routine, looping by chromosome
        for(String chromosome : oldVariants.getBreakendMap().keySet())
        {
            List<Breakend> oldBreakends = oldVariants.getChromosomeBreakends(chromosome);
            List<Breakend> newBreakends = newVariants.getChromosomeBreakends(chromosome);
            List<Breakend> oldFilteredBreakends = oldFilteredVariants.getChromosomeBreakends(chromosome);
            List<Breakend> newFilteredBreakends = newFilteredVariants.getChromosomeBreakends(chromosome);

            // first exact only from the main lists
            findBreakendMatches(oldBreakends, newBreakends, true);

            // then inexact from the main list
            findBreakendMatches(oldBreakends, newBreakends, false);

            findBreakendMatches(oldBreakends, newFilteredBreakends, false);

            findBreakendMatches(newBreakends, oldFilteredBreakends, false);

            // finally log unmatched breakends
            for(Breakend breakend : oldBreakends)
            {
                if(breakend.isFiltered() && !mConfig.IncludeNonPass)
                    continue;

                if(breakend.coordMatchType() == NONE)
                    writeVariantComparison(NONE, breakend.sv(), null, Collections.emptyList());
            }

            for(Breakend breakend : newBreakends)
            {
                if(breakend.isFiltered() && !mConfig.IncludeNonPass)
                    continue;

                if(breakend.coordMatchType() == NONE)
                    writeVariantComparison(NONE, null, breakend.sv(), Collections.emptyList());
            }
        }

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("VCF comparison complete");
    }

    private class BreakendMatch
    {
        public final Breakend Breakend;
        public final CoordMatchType MatchType;

        public BreakendMatch(final Breakend breakend, final CoordMatchType matchType)
        {
            Breakend = breakend;
            MatchType = matchType;
        }
    }

    private void findBreakendMatches(final List<Breakend> breakends, final List<Breakend> otherBreakends, boolean exactOnly)
    {
        for(Breakend breakend : breakends)
        {
            if(breakend.coordMatchType() != null)
                continue;

            if(breakend.isFiltered() && !mConfig.IncludeNonPass)
                continue;

            BreakendMatch match = findBreakendMatch(breakend, otherBreakends);

            if(match == null)
                continue;

            if(exactOnly && match.MatchType != EXACT)
                continue;

            breakend.setCoordMatchType(match.MatchType);
            match.Breakend.setCoordMatchType(match.MatchType);

            compareBreakends(breakend, match.Breakend, match.MatchType);
        }
    }

    private BreakendMatch findBreakendMatch(final Breakend breakend, final List<Breakend> otherBreakends)
    {
        CoordMatchType topMatchType = NONE;
        Breakend topBreakend = null;

        for(Breakend otherBreakend : otherBreakends)
        {
            if(otherBreakend.coordMatchType() != null)
                continue;

            CoordMatchType matchType = determineMatchType(breakend, otherBreakend);

            if(matchType.ordinal() < topMatchType.ordinal())
            {
                topMatchType = matchType;
                topBreakend = otherBreakend;
            }
        }

        return topBreakend != null ? new BreakendMatch(topBreakend, topMatchType) : null;
    }

    private void compareBreakends(final Breakend oldBreakend, final Breakend newBreakend, final CoordMatchType coordMatchType)
    {
        List<String> diffs = Lists.newArrayList();

        checkValueDifference(diffs, "Type", oldBreakend.type().toString(), newBreakend.type().toString());

        checkValueDifference(diffs, "Coords", oldBreakend.sv().svString(), newBreakend.sv().svString());

        checkValueDifference(
                diffs, "TumorFrags", oldBreakend.tumorFragmentCount(), newBreakend.tumorFragmentCount());

        checkValueDifference(
                diffs, "RefFrags", oldBreakend.referenceFragmentCount(), newBreakend.referenceFragmentCount());

        checkValueDifference(diffs, "Qual", oldBreakend.sv().qual(), newBreakend.sv().qual());

        checkValueDifference(diffs, "Passing", oldBreakend.isPass(), newBreakend.isPass());

        checkValueDifference(diffs, "InsSeqLength", oldBreakend.InsertSequence.length(), newBreakend.InsertSequence.length());

        checkValueDifference(diffs, "InexactHom", oldBreakend.InexactHomology.length(), newBreakend.InexactHomology.length());

        checkValueDifference(diffs, "Chained", oldBreakend.sv().inChainedAssembly(), newBreakend.sv().inChainedAssembly());

        checkValueDifference(diffs, "LINE", oldBreakend.sv().isLineSite(), newBreakend.sv().isLineSite());

        writeVariantComparison(coordMatchType, oldBreakend.sv(), newBreakend.sv(), diffs);
    }

    private BufferedWriter initialiseWriter(final String sampleId)
    {
        try
        {
            String fileName = mConfig.OutputDir + sampleId + "." + ESVEE_FILE_ID + ".vcf_compare";

            if(mConfig.OutputId != null)
                fileName += "." + mConfig.OutputId;

            fileName += TSV_EXTENSION;

            SV_LOGGER.debug("writing comparison file: {}", fileName);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(fileName, false);

            String header = String.join(
                    TSV_DELIM,
                    "SvCoords", "CoordMatch", "Diffs", "OldVcfId", "NewVcfId",
                    "OldTumorFrags",  "NewTumorFrags", "OldNormalFrags","NewNormalFrags", "OldQual", "NewQual", "OldFilters", "NewFilters");

            writer.write(header);
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            System.exit(1);
            return null;
        }
    }

    private void writeVariantComparison(
            final CoordMatchType coordMatchType, @Nullable final Variant oldVariant, @Nullable final Variant newVariant,
            final List<String> diffs)
    {
        if(diffs.isEmpty() && coordMatchType == EXACT && !mConfig.WriteMatches)
            return;

        if(oldVariant != null)
        {
            if(mComparedVariants.contains(oldVariant))
                return;

            mComparedVariants.add(oldVariant);
        }

        if(newVariant != null)
        {
            if(mComparedVariants.contains(newVariant))
                return;

            mComparedVariants.add(newVariant);
        }

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            String matchType = null;
            String varCoords = null;

            if(oldVariant != null && newVariant != null)
            {
                if(diffs.isEmpty())
                    matchType =  coordMatchType.toString();
                else
                    matchType = "VALUE";

                varCoords = oldVariant.svString();
            }
            else if(oldVariant != null)
            {
                matchType = "OLD_ONLY";
                varCoords = oldVariant.svString();
            }
            else
            {
                matchType = "NEW_ONLY";
                varCoords = newVariant.svString();
            }

            sj.add(varCoords).add(matchType);

            String diffStr = diffs.stream().collect(Collectors.joining(ITEM_DELIM));
            sj.add(diffStr);

            sj.add(oldVariant != null ? oldVariant.id() : "");
            sj.add(newVariant != null ? newVariant.id() : "");

            sj.add(String.valueOf(oldVariant != null ? oldVariant.breakendStart().tumorFragmentCount() : 0));
            sj.add(String.valueOf(newVariant != null ? newVariant.breakendStart().tumorFragmentCount() : 0));
            sj.add(String.valueOf(oldVariant != null ? oldVariant.breakendStart().referenceFragmentCount() : 0));
            sj.add(String.valueOf(newVariant != null ? newVariant.breakendStart().referenceFragmentCount() : 0));

            sj.add(String.valueOf(oldVariant != null ? oldVariant.qual() : 0));
            sj.add(String.valueOf(newVariant != null ? newVariant.qual() : 0));

            String oldFilters = oldVariant != null ? oldVariant.breakendStart().Filters.stream().collect(Collectors.joining(ITEM_DELIM)) : "";
            String newFilters = oldVariant != null ? newVariant.breakendStart().Filters.stream().collect(Collectors.joining(ITEM_DELIM)) : "";

            sj.add(oldFilters).add(newFilters);

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CompareConfig compareConfig = new CompareConfig(configBuilder);

        VcfComparer vcfComparer = new VcfComparer(compareConfig);
        vcfComparer.run();
    }
}
