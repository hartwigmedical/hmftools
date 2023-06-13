package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READ_COVERAGE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.depth.DepthConfig.VCF_TAG_PREFIX;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VcfDepthComparer
{
    private final String mVcfFile;
    private final String mNewVcfTagPrefix;
    private final double mDiffAbs;
    private final double mDiffPerc;

    private final BufferedWriter mWriter;

    private static final String VCF_FILE = "vcf_file";

    private static final int DEFAULT_MAX_DIFF = 2;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.02;

    public VcfDepthComparer(final CommandLine cmd)
    {
        mVcfFile = cmd.getOptionValue(VCF_FILE);
        mNewVcfTagPrefix = cmd.getOptionValue(VCF_TAG_PREFIX);
        mDiffAbs = DEFAULT_MAX_DIFF;
        mDiffPerc = DEFAULT_MAX_DIFF_PERC;

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mVcfFile == null || mNewVcfTagPrefix == null)
        {
            SV_LOGGER.error("missing VCF file or new tag prefix");
            return;
        }

        SV_LOGGER.info("loading VCF({})", mVcfFile);

        VcfFileReader reader = new VcfFileReader(mVcfFile);

        if(!reader.fileValid())
        {
            SV_LOGGER.error("error reading vcf({})", mVcfFile);
            System.exit(1);
        }

        List<String> oldVcfTags = Lists.newArrayList(REF_READ_COVERAGE, REF_READPAIR_COVERAGE);
        List<String> newVcfTags = oldVcfTags.stream().map(x -> format("%s_%s", mNewVcfTagPrefix, x)).collect(Collectors.toList());

        for(VariantContext variant : reader.iterator())
        {
            for(Genotype genotype : variant.getGenotypes())
            {
                if(genotype.getExtendedAttributes() == null || genotype.getExtendedAttributes().isEmpty())
                    continue;

                for(int i = 0; i < oldVcfTags.size(); ++i)
                {
                    String oldTag = oldVcfTags.get(i);
                    String newTag = newVcfTags.get(i);

                    int oldValue = getGenotypeAttributeAsInt(genotype, oldTag, 0);
                    int newValue = getGenotypeAttributeAsInt(genotype, newTag, 0);

                    if(hasDiff(oldValue, newValue))
                    {
                        writeDiffs(variant, oldTag, genotype.getSampleName(), oldValue, newValue);
                    }
                }
            }
        }

        closeBufferedWriter(mWriter);
    }

    private boolean hasDiff(int value1, int value2)
    {
        double diff = abs(value1 - value2);
        double diffPerc = diff / (double)max(value1, value2);
        return diff > mDiffAbs && diffPerc > mDiffPerc;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            int extensionIndex = mVcfFile.indexOf("vcf.gz");
            String fileName = mVcfFile.substring(0, extensionIndex) + "compare.csv";

            SV_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("VcfId,Chromosome,Position,SampleId,Tag,OrigValue,NewValue,Qual,Filters");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiffs(final VariantContext variant, final String vcfTag, final String sampleId, int origValue, int newValue)
    {
        try
        {
            mWriter.write(format("%s,%s,%d,%s,%s,%d,%d,%.1f,%s",
                    variant.getID(), variant.getContig(), variant.getStart(), sampleId, vcfTag, origValue, newValue,
                    variant.getPhredScaledQual(), variant.getFilters()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Optional, name of the reference sample");
        options.addOption(VCF_TAG_PREFIX, true, "New VCF tag prefix");

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        VcfDepthComparer vcfDepthComparer = new VcfDepthComparer(cmd);
        vcfDepthComparer.run();

        SV_LOGGER.info("RefDepth comparer complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
