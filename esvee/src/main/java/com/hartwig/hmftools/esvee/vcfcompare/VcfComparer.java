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
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFHeader;

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
        if(mConfig.SampleIds == null)
        {
            String sampleId = extractVcfSampleId(mConfig.OldVcf);
            mWriter = initialiseWriter(sampleId);

            SampleProcessor sampleProcessor = new SampleProcessor(sampleId);
            sampleProcessor.call();
        }
        else
        {
            mWriter = initialiseWriter(null);

            List<SampleProcessor> sampleTasks = Lists.newArrayList();

            for(String sampleId : mConfig.SampleIds)
            {
                SampleProcessor sampleProcessor = new SampleProcessor(sampleId);
                sampleTasks.add(sampleProcessor);
            }

            final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("VCF comparison complete");

    }

    private static String extractVcfSampleId(final String vcfFilename)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFilename);
        VCFHeader vcfHeader = vcfFileReader.vcfHeader();
        List<String> vcfSampleNames = vcfHeader.getGenotypeSamples();
        return vcfSampleNames.size() > 1 ? vcfSampleNames.get(1) :  vcfSampleNames.get(0);
    }

    private class SampleProcessor implements Callable<Void>
    {
        private final String mSampleId;

        public SampleProcessor(final String sampleId)
        {
            mSampleId = sampleId;
        }

        private static String formFilename(final String sampleId, final String filename)
        {
            return filename.replaceAll("\\*", sampleId);
        }

        @Override
        public Void call()
        {
           // SV_LOGGER.info("comparing VCFs for sample: " + mSampleId);
            VariantCache oldVariants = new VariantCache(formFilename(mSampleId, mConfig.OldVcf));
            VariantCache newVariants = new VariantCache(formFilename(mSampleId, mConfig.NewVcf));

            VariantCache oldFilteredVariants = new VariantCache(formFilename(mSampleId, mConfig.OldUnfilteredVcf));
            VariantCache newFilteredVariants = new VariantCache(formFilename(mSampleId, mConfig.NewUnfilteredVcf));

            Set<String> combinedChromosomes = Sets.newHashSet(oldVariants.getBreakendMap().keySet());
            combinedChromosomes.addAll(newVariants.getBreakendMap().keySet());

            // begin matching routine, looping by chromosome
            for(String chromosome : combinedChromosomes)
            {
                List<Breakend> oldBreakends = oldVariants.getChromosomeBreakends(chromosome);
                List<Breakend> newBreakends = newVariants.getChromosomeBreakends(chromosome);
                List<Breakend> oldFilteredBreakends = oldFilteredVariants.getChromosomeBreakends(chromosome);
                List<Breakend> newFilteredBreakends = newFilteredVariants.getChromosomeBreakends(chromosome);

                // first exact only from the passing breakeands
                findBreakendMatches(oldBreakends, newBreakends, true);

                // then inexact from the main lists
                findBreakendMatches(oldBreakends, newBreakends, false);

                // then search for old in the filtered new list and vice versa - again first exact then inexact matches
                findBreakendMatches(oldBreakends, newFilteredBreakends, true);
                findBreakendMatches(newBreakends, oldFilteredBreakends, true);
                findBreakendMatches(oldBreakends, newFilteredBreakends, false);
                findBreakendMatches(newBreakends, oldFilteredBreakends, false);

                // finally log unmatched breakends
                for(Breakend breakend : oldBreakends)
                {
                    if(breakend.isFiltered() && !mConfig.IncludeNonPass)
                        continue;

                    if(!breakend.matched())
                        writeVariantComparison(mSampleId, NONE, breakend.sv(), null, Collections.emptyList());
                }

                for(Breakend breakend : newBreakends)
                {
                    if(breakend.isFiltered() && !mConfig.IncludeNonPass)
                        continue;

                    if(!breakend.matched())
                        writeVariantComparison(mSampleId, NONE, null, breakend.sv(), Collections.emptyList());
                }
            }

            if(mConfig.isMultiSample())
                SV_LOGGER.debug("sample({}) comparison complete", mSampleId);

            return null;
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

        private void findBreakendMatches(
                final List<Breakend> breakends, final List<Breakend> otherBreakends, boolean exactOnly)
        {
            for(Breakend breakend : breakends)
            {
                if(breakend.matched())
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
                if(otherBreakend.matched())
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

            writeVariantComparison(mSampleId, coordMatchType, oldBreakend.sv(), newBreakend.sv(), diffs);
        }
    }

    private BufferedWriter initialiseWriter(@Nullable final String sampleId)
    {
        try
        {
            String fileName = mConfig.OutputDir;

            if(sampleId != null)
            {
                fileName += sampleId;
            }
            else
            {
                fileName = "cohort";
            }

            fileName += "." + ESVEE_FILE_ID + ".vcf_compare";

            if(mConfig.OutputId != null)
                fileName += "." + mConfig.OutputId;

            fileName += TSV_EXTENSION;

            SV_LOGGER.debug("writing comparison file: {}", fileName);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(sampleId == null)
            {
                sj.add("SampleId");
            }

            sj.add("SvType");
            sj.add("SvLength");

            sj.add("SvCoords");
            sj.add("CoordMatch");
            sj.add("Diffs");
            sj.add("OldVcfId");
            sj.add("NewVcfId");
            sj.add("OldTumorFrags");
            sj.add("NewTumorFrags");
            sj.add("OldNormalFrags");
            sj.add("NewNormalFrags");
            sj.add("OldQual");
            sj.add("NewQual");
            sj.add("OldFilters");
            sj.add("NewFilters");

            writer.write(sj.toString());
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

    private synchronized void writeVariantComparison(
            final String sampleId, final CoordMatchType coordMatchType, @Nullable final Variant oldVariant,
            @Nullable final Variant newVariant, final List<String> diffs)
    {
        if(diffs.isEmpty() && coordMatchType == EXACT && !mConfig.WriteMatches)
            return;

        // ignore non-passing matches
        if(oldVariant != null && newVariant != null && oldVariant.isFiltered() && newVariant.isFiltered())
        {
            return;
        }
        else if(oldVariant != null && oldVariant.isFiltered())
        {
            return;
        }
        else if(newVariant != null && newVariant.isFiltered())
        {
            return;
        }

        if(!mConfig.WritePonBreakends)
        {
            if((oldVariant == null || oldVariant.isPonOnly()) && (newVariant == null || newVariant.isPonOnly()))
                return;
        }

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

            if(mConfig.isMultiSample())
                sj.add(sampleId);

            String matchType, varCoords ;

            if(oldVariant != null && newVariant != null)
            {
                sj.add(String.valueOf(oldVariant.type()));
                sj.add(String.valueOf(oldVariant.length()));

                if(diffs.isEmpty())
                    matchType =  coordMatchType.toString();
                else
                    matchType = "VALUE";

                varCoords = oldVariant.svString();
            }
            else if(oldVariant != null)
            {
                sj.add(String.valueOf(oldVariant.type()));
                sj.add(String.valueOf(oldVariant.length()));

                matchType = "OLD_ONLY";
                varCoords = oldVariant.svString();
            }
            else
            {
                sj.add(String.valueOf(newVariant.type()));
                sj.add(String.valueOf(newVariant.length()));

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
            String newFilters = newVariant != null ? newVariant.breakendStart().Filters.stream().collect(Collectors.joining(ITEM_DELIM)) : "";

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
