package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.sage.SageApplication.createCommandLine;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.pipeline.ChromosomePartition;
import com.hartwig.hmftools.sage.quality.BaseQualityRecalibration;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class SageAppendApplication
{
    private final SageConfig mConfig;
    private final String mInputVcf;
    private final ExecutorService mExecutorService;
    private final IndexedFastaSequenceFile mRefGenome;

    private static final double MIN_PRIOR_VERSION = 3.0;
    private static final String INPUT_VCF = "input_vcf";

    public SageAppendApplication(final Options options, final String... args) throws ParseException, IOException
    {
        final VersionInfo version = new VersionInfo("sage.version");
        SG_LOGGER.info("SAGE version: {}", version.version());

        final CommandLine cmd = createCommandLine(args, options);
        mConfig = new SageConfig(true, version.version(), cmd);
        mInputVcf = mConfig.SampleDataDir + cmd.getOptionValue(INPUT_VCF);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
        mRefGenome = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));
    }

    public void run() throws IOException, ExecutionException, InterruptedException
    {
        // check config
        if(mInputVcf == null || mInputVcf.isEmpty())
        {
            SG_LOGGER.error("no input VCF file specified");
            System.exit(1);
        }

        if(mInputVcf.equals(mConfig.OutputFile))
        {
            SG_LOGGER.error("input and output VCFs must be different");
            System.exit(1);
        }

        if(mConfig.ReferenceIds.isEmpty())
        {
            SG_LOGGER.error("missing reference Id must be supplied");
            System.exit(1);
        }

        SG_LOGGER.info("reading and validating file: {}", mInputVcf);

        long startTime = System.currentTimeMillis();

        final AbstractFeatureReader<VariantContext, LineIterator> vcfReader = AbstractFeatureReader.getFeatureReader(
                mInputVcf, new VCFCodec(), false);

        VCFHeader inputHeader = (VCFHeader) vcfReader.getHeader();
        if(!validateInputHeader(inputHeader))
        {
            System.exit(1);
        }

        SG_LOGGER.info("writing to file: {}", mConfig.OutputFile);
        final VariantVCF outputVCF = new VariantVCF(mRefGenome, mConfig, inputHeader);

        final List<VariantContext> existingVariants = verifyAndReadExisting(vcfReader);

        if(existingVariants == null)
        {
            System.exit(1);
        }

        SG_LOGGER.info("loaded {} variants", existingVariants.size());

        final SAMSequenceDictionary dictionary = dictionary();
        // final List<Future<List<VariantContext>>> futures = Lists.newArrayList();

        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(mConfig, mExecutorService, mRefGenome);
        baseQualityRecalibration.produceRecalibrationMap();
        final Map<String,QualityRecalibrationMap> recalibrationMap = baseQualityRecalibration.getSampleRecalibrationMap();

        // final AdditionalReferencePipeline pipeline = new AdditionalReferencePipeline(mConfig, mExecutorService, mRefGenome, recalibrationMap);

        final ChromosomePartition chromosomePartition = new ChromosomePartition(mConfig, mRefGenome);

        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) && !MitochondrialChromosome.contains(chromosome))
                continue;

            SG_LOGGER.info("processing chromosome({})", chromosome);

            final List<VariantContext> chromosomeVariants = existingVariants.stream()
                    .filter(x -> x.getContig().equals(chromosome)).collect(Collectors.toList());

            List<ChrBaseRegion> chrBaseRegions = chromosomePartition.partition(chromosome);

            List<RegionAppendTask> regionTasks = Lists.newArrayList();

            for(int i = 0; i < chrBaseRegions.size(); ++i)
            {
                ChrBaseRegion region = chrBaseRegions.get(i);

                final List<VariantContext> regionVariants = chromosomeVariants.stream()
                        // .filter(x -> x.getStart() >= region.start() && x.getStart() <= region.end()) // would ignore variants spanning a region
                        .filter(x -> positionWithin(x.getStart(), region.start(), region.end()))
                        .collect(Collectors.toList());

                if(regionVariants.isEmpty())
                    continue;

                regionTasks.add(new RegionAppendTask(i, region, regionVariants, mConfig, mRefGenome, recalibrationMap));

                // futures.add(pipeline.appendReference(region, regionVariants));
            }

            final List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);

            for(RegionAppendTask regionTask : regionTasks)
            {
                final List<VariantContext> updatedVariants = regionTask.finalVariants();
                updatedVariants.forEach(outputVCF::write);
            }
        }

        /*
        for(Future<List<VariantContext>> updatedVariantsFuture : futures)
        {
            final List<VariantContext> updatedVariants = updatedVariantsFuture.get();
            updatedVariants.forEach(outputVCF::write);
        }
        */

        vcfReader.close();
        outputVCF.close();

        mRefGenome.close();

        long timeTaken = System.currentTimeMillis() - startTime;
        SG_LOGGER.info("completed in {} seconds", String.format("%.1f",timeTaken / 1000.0));
    }

    private List<VariantContext> verifyAndReadExisting(final AbstractFeatureReader<VariantContext, LineIterator> vcfReader)
    {
        try
        {
            List<VariantContext> result = Lists.newArrayList();

            VCFHeader header = (VCFHeader) vcfReader.getHeader();

            for(VariantContext variantContext : vcfReader.iterator())
            {
                result.add(variantContext.fullyDecode(header, false));
            }

            return result;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to read intput VCF: {}", e.toString());
            return null;
        }
    }

    private boolean validateInputHeader(VCFHeader header)
    {
        double oldVersion = sageVersion(header);
        if(Doubles.lessThan(oldVersion, MIN_PRIOR_VERSION))
        {
            SG_LOGGER.error("Sage VCF version({}) older than required({})", oldVersion, MIN_PRIOR_VERSION);
            return false;
        }

        final Set<String> existingSamples = existingSamples(header);

        StringJoiner sj = new StringJoiner(", ");
        existingSamples.forEach(x -> sj.add(x));

        SG_LOGGER.info("existing VCF samples: {}", sj.toString());

        for(String refSample : mConfig.ReferenceIds)
        {
            if(existingSamples.contains(refSample))
            {
                SG_LOGGER.error("config reference sample({}) already exits in input VCF", refSample);
                return false;
            }
        }

        return true;
    }

    private static double sageVersion(@NotNull final VCFHeader header)
    {
        VCFHeaderLine oldVersion = header.getMetaDataLine(VariantVCF.VERSION_META_DATA);

        if(oldVersion == null)
            return 0;

        String oldVersionString = oldVersion.getValue();
        try
        {
            return Double.parseDouble(oldVersionString);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to parse Sage version: {}", oldVersionString);
            return 0;
        }
    }

    private static Set<String> existingSamples(final VCFHeader header)
    {
        return Sets.newHashSet(header.getGenotypeSamples());
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = mConfig.ReferenceBams.get(0);

        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome)).open(new File(bam));

        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

    public static Options createOptions()
    {
        final Options options = new Options();
        SageConfig.commonOptions().getOptions().forEach(options::addOption);
        options.addOption(INPUT_VCF, true, "Path to input vcf");
        return options;
    }

    public static void main(String[] args)
    {
        final Options options = createOptions();

        try
        {
            final SageAppendApplication application = new SageAppendApplication(options, args);
            application.run();
        }
        catch(ParseException e)
        {
            SG_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageAppendApplication", options);
            System.exit(1);
        }
        catch(Exception e)
        {
            SG_LOGGER.warn(e);
            System.exit(1);
        }
    }

}
