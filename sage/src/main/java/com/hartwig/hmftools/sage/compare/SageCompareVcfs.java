package com.hartwig.hmftools.sage.compare;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Set;
import java.util.StringJoiner;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

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

    private final BufferedWriter mWriter;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";
    private static final String SAMPLE = "sample";

    public SageCompareVcfs(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mOriginalVcf = cmd.getOptionValue(ORIGINAL_VCF);
        mNewVcf = cmd.getOptionValue(NEW_VCF);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SG_LOGGER.error("missing VCFs");
            return;
        }

        SG_LOGGER.info("comparing VCFs: orig({}) vs new({})", mOriginalVcf, mNewVcf);

        AbstractFeatureReader<VariantContext, LineIterator> origReader = AbstractFeatureReader.getFeatureReader(mOriginalVcf, new VCFCodec(), false);
        AbstractFeatureReader<VariantContext, LineIterator> newReader = AbstractFeatureReader.getFeatureReader(mNewVcf, new VCFCodec(), false);

        int unmatchedOrigCount = 0;
        int unmatchedNewCount = 0;
        int compareCount = 0;

        try
        {
            Iterator origIter = origReader.iterator();
            Iterator newIter = newReader.iterator();

            VariantCompareData origVar = VariantCompareData.fromContext((VariantContext)origIter.next());
            VariantCompareData newVar = VariantCompareData.fromContext((VariantContext)newIter.next());

            while(newVar != null || origVar != null)
            {
                int totalComparisons = compareCount + unmatchedNewCount + unmatchedOrigCount;

                if(totalComparisons > 0 && (totalComparisons % 1000) == 0)
                {
                    SG_LOGGER.info("processed {} variants", totalComparisons);
                }

                if(newVar != null && origVar != null)
                {
                    int posCompare = comparePositions(origVar, newVar);

                    if(posCompare == 0)
                    {
                        compareVariants(origVar, newVar);
                        ++compareCount;

                        origVar = getNextVariant(origIter);
                        newVar = getNextVariant(newIter);
                    }
                    else if(posCompare < 0)
                    {
                        writeDiffs(origVar, null, "NO_NEW", "", "");
                        ++unmatchedOrigCount;

                        origVar = getNextVariant(origIter);
                    }
                    else
                    {
                        writeDiffs(null, newVar, "NO_ORIG", "", "");
                        ++unmatchedNewCount;

                        newVar = getNextVariant(newIter);
                    }
                }
                else if(newVar != null)
                {
                    writeDiffs(null, newVar, "NO_ORIG", "", "");
                    ++unmatchedNewCount;

                    newVar = getNextVariant(newIter);
                }
                else if(origVar != null)
                {
                    writeDiffs(origVar, null, "NO_NEW", "", "");
                    ++unmatchedOrigCount;

                    origVar = getNextVariant(origIter);
                }
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("error reading vcf files: {}", e.toString());
        }

        SG_LOGGER.info("summary: compared({}) unmatched(orig={} new={})", compareCount, unmatchedOrigCount, unmatchedNewCount);

        closeBufferedWriter(mWriter);
    }

    private static VariantCompareData getNextVariant(final Iterator iter)
    {
        // iter.next();
        return VariantCompareData.fromContext((VariantContext)iter.next());
    }

    private static int comparePositions(final VariantCompareData first, final VariantCompareData second)
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

    private void compareVariants(final VariantCompareData origVar, final VariantCompareData newVar)
    {
        // compare quals
        // writeDiffs


        // compare filters



        // compare depth and any other critical info

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

            writer.write("Chromosome,Position,Ref,Alt,Tier,DiffType,OrigValue,NewValue,OrigQual,NewQual,OrigFilters,NewFilters");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiffs(
            final VariantCompareData origVar, final VariantCompareData newVar,
            final String diffType, final String origValue, final String newValue)
    {
        try
        {
            VariantCompareData var = origVar != null ? origVar : newVar;

            mWriter.write(String.format("%s,%d,%s,%s,%s",
                    var.Chromosome, var.Position, var.Ref, var.Alt, var.tier()));

            mWriter.write(String.format(",%s,%s,%s",
                    diffType, origValue, newValue));

            mWriter.write(String.format(",%.1f,%.1f,%s,%s",
                    origVar != null ? origVar.qual() : -1,
                    newVar != null ? newVar.qual() : -1,
                    origVar != null ? filtersStr(origVar.filters()) : "",
                    newVar != null ? filtersStr(newVar.filters()) : ""));

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
