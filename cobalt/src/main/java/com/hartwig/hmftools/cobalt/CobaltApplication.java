package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;
import com.google.common.collect.Multimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.count.CountSupplier;
import com.hartwig.hmftools.cobalt.targeted.TargetRegionEnrichment;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class CobaltApplication implements AutoCloseable
{
    @ParametersDelegate
    private final CobaltConfig mConfig = new CobaltConfig();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        CB_LOGGER.info("{}", LocalDate.now());
        CB_LOGGER.info("args: {}", String.join(" ", args));

        final CobaltApplication application = new CobaltApplication();
        JCommander commander = JCommander.newBuilder()
                .addObject(application)
                .build();

        // use unix style formatter
        commander.setUsageFormatter(new UnixStyleUsageFormatter(commander));
        // help message show in order parameters are declared
        commander.setParameterDescriptionComparator(new DeclaredOrderParameterComparator(application.getClass()));

        try
        {
            commander.parse(args);
            System.exit(application.run());
        }
        catch (com.beust.jcommander.ParameterException e)
        {
            System.out.println("Unable to parse args: " + e.getMessage());
            commander.usage();
            System.exit(1);
        }
    }

    private int run() throws IOException, ExecutionException, InterruptedException
    {
        mConfig.validate();
        mLoggingOptions.setLogLevel();

        VersionInfo mVersionInfo = new VersionInfo("cobalt.version");
        CB_LOGGER.info("COBALT version: {}", mVersionInfo.version());

        CB_LOGGER.info("Reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        final SamReaderFactory readerFactory = readerFactory(mConfig);

        final Collection<Chromosome> chromosomes = loadChromosomes(readerFactory);

        final Multimap<Chromosome, GCProfile> gcProfiles = loadGCContent(chromosomes);

        final CountSupplier countSupplier = new CountSupplier(
                WINDOW_SIZE, mConfig.MinMappingQuality,
                executorService, readerFactory, chromosomes);
        
        countSupplier.generateCounts(mConfig.ReferenceBamPath, mConfig.TumorBamPath);

        final RatioSupplier ratioSupplier = new RatioSupplier(mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir,
                gcProfiles, chromosomes, countSupplier.getReferenceCounts(), countSupplier.getTumorCounts());

        if (mConfig.TargetRegionPath != null)
        {
            var targetRegionEnrichment = TargetRegionEnrichment.fromTsv(mConfig.TargetRegionPath);
            ratioSupplier.setTargetRegionEnrichment(targetRegionEnrichment);
        }

        Multimap<Chromosome, CobaltRatio> ratios;

        switch (mConfig.mode())
        {
            case TUMOR_ONLY:
                ratios = ratioSupplier.tumorOnly(mConfig.TumorOnlyDiploidBed);
                break;
            case GERMLIHE_ONLY:
                ratios = ratioSupplier.germlineOnly();
                break;
            default:
                ratios = ratioSupplier.tumorNormalPair();
        }

        final String outputFilename = CobaltRatioFile.generateFilenameForWriting(
                mConfig.OutputDir, mConfig.TumorId != null ? mConfig.TumorId : mConfig.ReferenceId);
        CB_LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
        mVersionInfo.write(mConfig.OutputDir);
        CobaltRatioFile.write(outputFilename, ratios.values());

        applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);

        executorService.shutdown();

        CB_LOGGER.info("Complete");

        return 0;
    }

    @NotNull
    private static SamReaderFactory readerFactory(@NotNull final CobaltConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.Stringency);
        if(config.RefGenomePath != null && !config.RefGenomePath.isEmpty())
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomePath)));
        }
        return readerFactory;
    }

    @Override
    public void close()
    {
        CB_LOGGER.info("Complete");
    }

    private Collection<Chromosome> loadChromosomes(final SamReaderFactory readerFactory) throws IOException
    {
        Collection<Chromosome> chromosomes = new ArrayList<>();

        try(SamReader reader = readerFactory.open(new File(mConfig.ReferenceBamPath != null ? mConfig.ReferenceBamPath : mConfig.TumorBamPath)))
        {
            final List<ChromosomeLength> chromosomeLengths = ChromosomeLengthFactory.create(reader.getFileHeader());
            chromosomeLengths.forEach(o -> chromosomes.add(new Chromosome(o.chromosome(), o.length())));
        }

        return chromosomes;
    }

    @NotNull
    public Multimap<Chromosome, GCProfile> loadGCContent(Collection<Chromosome> chromosomes) throws IOException
    {
        Collection<GCProfile> gcProfileList = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath).values();
        return CobaltUtils.toChromosomeMultiMap(gcProfileList, chromosomes, GCProfile::chromosome);
    }
}
