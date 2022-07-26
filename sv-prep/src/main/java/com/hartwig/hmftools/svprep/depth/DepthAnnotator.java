package com.hartwig.hmftools.svprep.depth;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.SAMPLE;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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

public class DepthAnnotator
{
    private final String mInputVcf;
    private final String mOutputVcf;
    private final String mBamFile;
    private final String mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final List<ChrBaseRegion> mSpecificRegions;

    private final int mThreads;

    private final Map<String, List<VariantContext>> mChrVariantMap;

    private static final String INPUT_VCF = "input_vcf";
    private static final String OUTPUT_VCF = "output_vcf";
    private static final String BAM_FILE = "bam_file";
    private static final String THREADS = "threads";

    public DepthAnnotator(final CommandLine cmd)
    {
        mInputVcf = cmd.getOptionValue(INPUT_VCF);
        mOutputVcf = cmd.getOptionValue(OUTPUT_VCF);
        mBamFile = cmd.getOptionValue(BAM_FILE);
        mRefGenome = cmd.getOptionValue(REF_GENOME);
        mRefGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));;
        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        mSpecificRegions = Lists.newArrayList();

        try
        {
            mSpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }

        mChrVariantMap = Maps.newHashMap();
    }

    public void run()
    {
        if(mInputVcf == null || mOutputVcf == null)
        {
            SV_LOGGER.error("missing VCF");
            System.exit(1);
        }

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                mInputVcf, new VCFCodec(), false);

        int vcfCount = 0;

        try
        {
            String currentChromosome = "";
            List<VariantContext> variantsList = null;

            for(VariantContext variantContext : reader.iterator())
            {
                ++vcfCount;

                if(excludeVariant(variantContext))
                    continue;

                String chromosome = variantContext.getContig();

                if(!chromosome.equals(currentChromosome))
                {
                    variantsList = Lists.newArrayList();
                    mChrVariantMap.put(chromosome, variantsList);
                    currentChromosome = chromosome;
                }

                variantsList.add(variantContext);
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", mInputVcf, e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("loaded {} variants from vcf({})", vcfCount, mInputVcf);

        List<DepthTask> depthTasks = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            List<VariantContext> variantsList = mChrVariantMap.get(chrStr);

            if(variantsList == null)
                continue;

            DepthTask depthTask = new DepthTask(chromosome.toString(), mRefGenome, mBamFile);
            depthTask.variants().addAll(variantsList);
            depthTasks.add(depthTask);
        }

        final List<Callable> callableList = depthTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        SV_LOGGER.info("depth annotation complete");

        // write output VCF
    }

    private boolean excludeVariant(final VariantContext variant)
    {
        if(mSpecificRegions.isEmpty())
            return false;

        return mSpecificRegions.stream().noneMatch(x -> x.containsPosition(variant.getContig(), variant.getStart()));
    }


    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the sample");
        options.addOption(INPUT_VCF, true, "Input VCF File");
        options.addOption(OUTPUT_VCF, true, "Output VCF File");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(BAM_FILE, true, "BAM file to slice for depth");
        options.addOption(THREADS, true, "Multi-thread count");

        addSpecificChromosomesRegionsConfig(options);
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        DepthAnnotator bamSvSlicer = new DepthAnnotator(cmd);
        bamSvSlicer.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
