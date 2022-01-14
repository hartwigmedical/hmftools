package com.hartwig.hmftools.sage.compare;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sage.SageMetaData.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.compare.DiffType.ALLELE_DEPTH;
import static com.hartwig.hmftools.sage.compare.DiffType.FILTER_DIFF;
import static com.hartwig.hmftools.sage.compare.DiffType.FILTER_PASS;
import static com.hartwig.hmftools.sage.compare.DiffType.LOCAL_PHASE;
import static com.hartwig.hmftools.sage.compare.DiffType.NO_NEW;
import static com.hartwig.hmftools.sage.compare.DiffType.NO_ORIG;
import static com.hartwig.hmftools.sage.compare.DiffType.OTHER_VALUE;
import static com.hartwig.hmftools.sage.compare.DiffType.QUAL;
import static com.hartwig.hmftools.sage.compare.DiffType.TIER;
import static com.hartwig.hmftools.sage.compare.DiffType.hasValueDiff;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.PASS;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

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

    private final BufferedWriter mWriter;

    // counters and state
    private int mCompareCount;
    private int mDiffCount;
    private int mUnmatchedOrigCount;
    private int mUnmatchedNewCount;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";
    private static final String SAMPLE = "sample";

    private static final String DIFF_ABS = "diff_abs";
    private static final String DIFF_PERC = "diff_perc";

    private static final double DEFAULT_DIFF_PERC = 0.1;
    private static final double DEFAULT_DIFF_ABS = 2;

    public SageCompareVcfs(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mOriginalVcf = cmd.getOptionValue(ORIGINAL_VCF);
        mNewVcf = cmd.getOptionValue(NEW_VCF);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mDiffAbsPermitted = Double.parseDouble(cmd.getOptionValue(DIFF_ABS, String.valueOf(DEFAULT_DIFF_ABS)));
        mDiffPercPermitted = Double.parseDouble(cmd.getOptionValue(DIFF_PERC, String.valueOf(DEFAULT_DIFF_PERC)));

        mWriter = initialiseWriter();

        mDiffCount = 0;
        mCompareCount = 0;
        mUnmatchedOrigCount = 0;
        mUnmatchedNewCount = 0;
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SG_LOGGER.error("missing VCFs in config");
            System.exit(1);
        }

        if(!Files.exists(Paths.get(mOriginalVcf)) || !Files.exists(Paths.get(mNewVcf)))
        {
            SG_LOGGER.error("invalid VCF paths: original({}) and new({})", mOriginalVcf, mNewVcf);
            System.exit(1);
        }

        SG_LOGGER.info("comparing VCFs: orig({}) vs new({})", mOriginalVcf, mNewVcf);

        AbstractFeatureReader<VariantContext, LineIterator> origReader = AbstractFeatureReader.getFeatureReader(mOriginalVcf, new VCFCodec(), false);
        AbstractFeatureReader<VariantContext, LineIterator> newReader = AbstractFeatureReader.getFeatureReader(mNewVcf, new VCFCodec(), false);

        try
        {
            Iterator origIter = origReader.iterator();
            Iterator newIter = newReader.iterator();

            VariantData origVar = VariantData.fromContext((VariantContext)origIter.next());
            VariantData newVar = VariantData.fromContext((VariantContext)newIter.next());

            while(newVar != null || origVar != null)
            {
                int totalComparisons = totalComparisons();

                if(totalComparisons > 0 && (totalComparisons % 10000) == 0)
                {
                    SG_LOGGER.info("processed {} variants", totalComparisons);
                }

                if(newVar != null && origVar != null)
                {
                    int posCompare = comparePositions(origVar, newVar);

                    if(posCompare == 0)
                    {
                        compareVariants(origVar, newVar);

                        origVar = getNextVariant(origIter);
                        newVar = getNextVariant(newIter);
                    }
                    else if(posCompare < 0)
                    {
                        writeUnmatchedVariant(newVar, false);
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
                    writeUnmatchedVariant(newVar, false);
                    origVar = getNextVariant(origIter);
                }
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("error reading vcf files: {}", e.toString());
        }

        SG_LOGGER.info("summary: total variants({}) compared({} diffs={}) unmatched(orig={} new={})",
                totalComparisons(), mCompareCount, mDiffCount, mUnmatchedOrigCount, mUnmatchedNewCount);

        closeBufferedWriter(mWriter);
    }

    private static VariantData getNextVariant(final Iterator iter)
    {
        return VariantData.fromContext((VariantContext)iter.next());
    }

    private int totalComparisons() { return mCompareCount + mUnmatchedNewCount + mUnmatchedOrigCount; }

    private static int comparePositions(final VariantData first, final VariantData second)
    {
        // return -1 if first is earlier in genomic space than second
        if(first.Chromosome.equals(second.Chromosome))
        {
            if(first.Position == second.Position)
                return 0;

            return first.Position < second.Position ? -1: 1;
        }
        else
        {
            return HumanChromosome.lowerChromosome(first.Chromosome, second.Chromosome) ? -1 : 1;
        }
    }

    private void compareVariants(final VariantData origVar, final VariantData newVar)
    {
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
            boolean origIsPass = origVar.context().isNotFiltered() || (origFilters.size() == 1 && origFilters.contains(PASS));
            boolean newIsPass = newVar.context().isNotFiltered() || (newFilters.size() == 1 && newFilters.contains(PASS));

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

            fileName += ".csv";

            SG_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Chromosome,Position,Ref,Alt,Tier,DiffType,OrigValue,NewValue,OrigQual,NewQual,SharedFilters,MaxDepth");
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
        if(isNew)
        {
            mUnmatchedNewCount++;
            writeDiffs(null, var, NO_ORIG, "", "");
        }
        else
        {
            mUnmatchedOrigCount++;
            writeDiffs(var, var, NO_NEW, "", "");
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

            mWriter.write(String.format("%s,%d,%s,%s,%s", var.Chromosome, var.Position, var.Ref, var.Alt, var.tier()));

            mWriter.write(String.format(",%s,%s,%s", diffType, origValue, newValue));

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

            mWriter.write(String.format(",%.0f,%.0f,%s,%s,%d",
                    origVar != null ? origVar.qual() : -1, newVar != null ? newVar.qual() : -1, sharedFilters, maxDepth));

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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(ORIGINAL_VCF, true, "Optional, name of the reference sample");
        options.addOption(NEW_VCF, true, "Path to the GRIDSS structural variant VCF file");

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        SageCompareVcfs sageCompare = new SageCompareVcfs(cmd);
        sageCompare.run();

        SG_LOGGER.info("Sage compare VCFs complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
