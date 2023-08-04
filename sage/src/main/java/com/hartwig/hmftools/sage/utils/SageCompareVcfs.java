package com.hartwig.hmftools.sage.utils;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.utils.DiffType.ALLELE_DEPTH;
import static com.hartwig.hmftools.sage.utils.DiffType.FILTER_DIFF;
import static com.hartwig.hmftools.sage.utils.DiffType.FILTER_PASS;
import static com.hartwig.hmftools.sage.utils.DiffType.LOCAL_PHASE;
import static com.hartwig.hmftools.sage.utils.DiffType.NO_NEW;
import static com.hartwig.hmftools.sage.utils.DiffType.NO_ORIG;
import static com.hartwig.hmftools.sage.utils.DiffType.QUAL;
import static com.hartwig.hmftools.sage.utils.DiffType.TIER;
import static com.hartwig.hmftools.sage.utils.DiffType.hasValueDiff;
import static com.hartwig.hmftools.sage.utils.VariantData.comparePositions;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SageCompareVcfs
{
    private final String mSampleId;
    private final String mOriginalVcf;
    private final String mNewVcf;
    private final String mOutputDir;
    private final String mOutputId;

    // comparison config
    private final double mDiffAbsPermitted;
    private final double mDiffPercPermitted;
    private final boolean mPassOnly;

    private final BufferedWriter mWriter;

    // counters and state
    private int mCompareCount;
    private int mDiffCount;
    private int mUnmatchedOrigCount;
    private int mUnmatchedNewCount;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";

    private static final String DIFF_ABS = "diff_abs";
    private static final String DIFF_PERC = "diff_perc";
    private static final String PASS_ONLY = "pass_only";

    private static final double DEFAULT_DIFF_PERC = 0.1;
    private static final double DEFAULT_DIFF_ABS = 2;

    public SageCompareVcfs(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mOriginalVcf = configBuilder.getValue(ORIGINAL_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mDiffAbsPermitted = Double.parseDouble(configBuilder.getValue(DIFF_ABS, String.valueOf(DEFAULT_DIFF_ABS)));
        mDiffPercPermitted = Double.parseDouble(configBuilder.getValue(DIFF_PERC, String.valueOf(DEFAULT_DIFF_PERC)));
        mPassOnly = configBuilder.hasFlag(PASS_ONLY);

        mWriter = initialiseWriter();

        mDiffCount = 0;
        mCompareCount = 0;
        mUnmatchedOrigCount = 0;
        mUnmatchedNewCount = 0;
    }

    private static final int LOG_COUNT = 100000;

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SG_LOGGER.error("missing VCFs in config");
            System.exit(1);
        }

        VcfFileReader origReader = new VcfFileReader(mOriginalVcf);
        VcfFileReader newReader = new VcfFileReader(mNewVcf);

        if(!origReader.fileValid() || !newReader.fileValid())
        {
            SG_LOGGER.error("invalid VCF paths: original({}) and new({})", mOriginalVcf, mNewVcf);
            System.exit(1);
        }

        SG_LOGGER.info("comparing VCFs: orig({}) vs new({})", mOriginalVcf, mNewVcf);

        Iterator origIter = origReader.iterator();
        Iterator newIter = newReader.iterator();

        VariantData origVar = VariantData.fromContext((VariantContext)origIter.next());
        VariantData newVar = VariantData.fromContext((VariantContext)newIter.next());

        List<VariantData> origVariants = Lists.newArrayList();
        List<VariantData> newVariants = Lists.newArrayList();

        int nextLog = LOG_COUNT;

        while(newVar != null || origVar != null)
        {
            int totalComparisons = totalComparisons();

            if(totalComparisons >= nextLog)
            {
                nextLog += LOG_COUNT;
                SG_LOGGER.info("processed {} variants", totalComparisons);
            }

            if(!origVariants.isEmpty() && !newVariants.isEmpty())
            {
                boolean positionMatch = false;

                if(origVar != null && comparePositions(origVariants.get(0), origVar) == 0)
                {
                    origVariants.add(origVar);
                    origVar = getNextVariant(origIter);
                    positionMatch = true;
                }

                if(newVar != null && comparePositions(newVariants.get(0), newVar) == 0)
                {
                    newVariants.add(newVar);
                    positionMatch = true;
                    newVar = getNextVariant(newIter);
                }

                if(positionMatch)
                    continue; // keep looking for more at this position

                // run comparisons within this group
                compareVariants(origVariants, newVariants);
                continue;
            }

            if(newVar != null && origVar != null)
            {
                int posCompare = comparePositions(origVar, newVar);

                if(posCompare == 0)
                {
                    origVariants.add(origVar);
                    newVariants.add(newVar);

                    // check for others at this exact position

                    origVar = getNextVariant(origIter);
                    newVar = getNextVariant(newIter);
                }
                else if(posCompare < 0)
                {
                    // original SV has a lower positions
                    writeUnmatchedVariant(origVar, false);
                    origVar = getNextVariant(origIter);
                }
                else
                {
                    writeUnmatchedVariant(newVar, true);
                    newVar = getNextVariant(newIter);
                }
            }
            else if(newVar != null)
            {
                writeUnmatchedVariant(newVar, true);
                newVar = getNextVariant(newIter);
            }
            else if(origVar != null)
            {
                writeUnmatchedVariant(origVar, false);
                origVar = getNextVariant(origIter);
            }
        }

        SG_LOGGER.info("summary: total variants({}) compared({} diffs={}) unmatched(orig={} new={})",
                totalComparisons(), mCompareCount, mDiffCount, mUnmatchedOrigCount, mUnmatchedNewCount);

        closeBufferedWriter(mWriter);

        SG_LOGGER.info("Sage compare VCFs complete");
    }

    private static VariantData getNextVariant(final Iterator iter)
    {
        return VariantData.fromContext((VariantContext)iter.next());
    }

    private int totalComparisons() { return mCompareCount + mUnmatchedNewCount + mUnmatchedOrigCount; }

    private void compareVariants(final List<VariantData> origVariants, final List<VariantData> newVariants)
    {
        int origIndex = 0;
        while(origIndex < origVariants.size())
        {
            VariantData origVar = origVariants.get(origIndex);

            boolean matched = false;
            int newIndex = 0;
            while(newIndex < newVariants.size())
            {
                VariantData newVar = newVariants.get(newIndex);

                if(origVar.matches(newVar))
                {
                    compareVariants(origVar, newVar);
                    newVariants.remove(newIndex);
                    origVariants.remove(origIndex);
                    matched = true;
                    break;
                }

                ++newIndex;
            }

            if(!matched)
                ++origIndex;
        }

        origVariants.forEach(x -> writeUnmatchedVariant(x, false));
        newVariants.forEach(x -> writeUnmatchedVariant(x, true));

        origVariants.clear();
        newVariants.clear();
    }

    private void compareVariants(final VariantData origVar, final VariantData newVar)
    {
        if(mPassOnly)
        {
            if(!origVar.isPassing() && !newVar.isPassing())
                return;

            if(!origVar.isPassing())
            {
                writeUnmatchedVariant(newVar, true);
                return;
            }
            else if(!newVar.isPassing())
            {
                writeUnmatchedVariant(origVar, false);
                return;
            }
        }

        mCompareCount++;

        // compare quals
        if(hasValueDiff(origVar.qual(), newVar.qual(), mDiffAbsPermitted, mDiffPercPermitted))
        {
            writeDiffs(origVar, newVar, QUAL, String.valueOf(origVar.qual()), String.valueOf(newVar.qual()));
            return;
        }

        if(origVar.tier() != newVar.tier())
        {
            writeDiffs(origVar, newVar, TIER, origVar.tier().toString(), newVar.tier().toString());
            return;
        }

        // compare filters

        Set<String> origFilters = origVar.filters();
        Set<String> newFilters = newVar.filters();

        // first check for a difference in PASS vs not
        Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
        Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

        if(!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty())
        {
            boolean origIsPass = origVar.isPassing();
            boolean newIsPass = newVar.isPassing();

            if(origIsPass != newIsPass)
            {
                writeDiffs(origVar, newVar, FILTER_PASS, filtersStr(origFilterDiffs), filtersStr(newFilterDiffs));
            }
            else
            {
                writeDiffs(origVar, newVar, FILTER_DIFF, filtersStr(origFilterDiffs), filtersStr(newFilterDiffs));
            }

            return;
        }

        // compare allelic depth
        int origAd = origVar.allelicDepth();
        int newAd = newVar.allelicDepth();

        if(hasValueDiff(origAd, newAd, mDiffAbsPermitted, mDiffPercPermitted))
        {
            writeDiffs(origVar, newVar, ALLELE_DEPTH, String.valueOf(origAd), String.valueOf(newAd));
            return;
        }

        boolean origPhased = origVar.context().hasAttribute(LOCAL_PHASE_SET);
        boolean newPhased = origVar.context().hasAttribute(LOCAL_PHASE_SET);

        if(origPhased != newPhased)
        {
            writeDiffs(origVar, newVar, LOCAL_PHASE, String.valueOf(origPhased), String.valueOf(newPhased));
            return;
        }

        // any other critical info
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".compare";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += TSV_EXTENSION;

            SG_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Chromosome\tPosition\tRef\tAlt\tTier\tDiffType\tOrigValue\tNewValue\tOrigQual\tNewQual\tSharedFilters\tMaxDepth");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeUnmatchedVariant(final VariantData var, boolean isNew)
    {
        if(mPassOnly && !var.isPassing())
            return;

        if(isNew)
        {
            mUnmatchedNewCount++;
            writeDiffs(null, var, NO_ORIG, "", "");
        }
        else
        {
            mUnmatchedOrigCount++;
            writeDiffs(var, null, NO_NEW, "", "");
        }
    }

    private void writeDiffs(
            final VariantData origVar, final VariantData newVar,
            final DiffType diffType, final String origValue, final String newValue)
    {
        if(origVar != null && newVar != null)
        {
            ++mDiffCount;
        }

        try
        {
            VariantData var = origVar != null ? origVar : newVar;

            mWriter.write(String.format("%s\t%d\t%s\t%s\t%s", var.Chromosome, var.Position, var.Ref, var.Alt, var.tier()));

            mWriter.write(String.format("\t%s\t%s\t%s", diffType, origValue, newValue));

            Set<String> sharedFilters = Sets.newHashSet();
            int maxDepth = 0;

            if(origVar != null)
            {
                origVar.filters().forEach(x -> sharedFilters.add(x));
                maxDepth = origVar.allelicDepth() + origVar.referenceDepth();
            }

            if(newVar != null)
            {
                newVar.filters().forEach(x -> sharedFilters.add(x));
                maxDepth = max(maxDepth, newVar.allelicDepth() + newVar.referenceDepth());
            }

            if(sharedFilters.isEmpty())
                sharedFilters.add(PASS);

            mWriter.write(String.format("\t%.0f\t%.0f\t%s\t%d",
                    origVar != null ? origVar.qual() : -1, newVar != null ? newVar.qual() : -1, filtersStr(sharedFilters), maxDepth));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private String filtersStr(final Set<String> filters)
    {
        StringJoiner sj = new StringJoiner(";");
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(ORIGINAL_VCF, true, "Original VCF file");
        configBuilder.addPath(NEW_VCF, true, "New VCF file");
        configBuilder.addFlag(PASS_ONLY, "Only compare passing variants");
        configBuilder.addDecimal(DIFF_ABS, "Absolute value difference", DEFAULT_DIFF_ABS);
        configBuilder.addDecimal(DIFF_PERC, "Percentage value difference", DEFAULT_DIFF_PERC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SageCompareVcfs sageCompare = new SageCompareVcfs(configBuilder);
        sageCompare.run();
    }
}
