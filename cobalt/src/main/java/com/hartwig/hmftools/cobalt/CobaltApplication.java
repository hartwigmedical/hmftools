package com.hartwig.hmftools.cobalt;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.cobalt.CobaltUtils.rowToCobaltRatio;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.Collection;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.count.BamReadCounter;
import com.hartwig.hmftools.cobalt.diploid.DiploidRegionLoader;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import tech.tablesaw.api.*;
import tech.tablesaw.io.csv.CsvReadOptions;

public class CobaltApplication
{
    @ParametersDelegate
    private final CobaltConfig mConfig = new CobaltConfig();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        // CB_LOGGER.debug("args: {}", String.join(" ", args));

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
        }
        catch (com.beust.jcommander.ParameterException e)
        {
            System.out.println("Unable to parse args: " + e.getMessage());
            commander.usage();
            System.exit(1);
        }

        // set all thread exception handler
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            CB_LOGGER.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(application.run(args));
    }

    private int run(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        Instant start = Instant.now();
        mConfig.validate();
        mLoggingOptions.setLogLevel();

        VersionInfo mVersionInfo = new VersionInfo("cobalt.version");

        CB_LOGGER.info("Cobalt version: {}", mVersionInfo.version());

        CB_LOGGER.debug("build timestamp: {}, run args: {}",
                mVersionInfo.buildTime().format(ISO_ZONED_DATE_TIME), String.join(" ", args));

        CB_LOGGER.info("Reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        try
        {
            final SamReaderFactory readerFactory = readerFactory(mConfig);

            var chromosomePosCodec = new ChromosomePositionCodec();

            final BamReadCounter bamReadCounter = new BamReadCounter(
                    WINDOW_SIZE, mConfig.MinMappingQuality,
                    executorService, readerFactory, chromosomePosCodec);

            bamReadCounter.generateDepths(mConfig.ReferenceBamPath, mConfig.TumorBamPath);

            Table referenceReadDepths = bamReadCounter.getReferenceDepths();
            Table tumorReadDepths = bamReadCounter.getTumorDepths();

            final Table gcProfiles = loadGCContent(chromosomePosCodec);

            final RatioSupplier ratioSupplier = new RatioSupplier(mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir,
                    gcProfiles, referenceReadDepths, tumorReadDepths,
                    chromosomePosCodec);

            if (mConfig.TargetRegionPath != null)
            {
                CsvReadOptions options = CsvReadOptions.builder(mConfig.TargetRegionPath).separator('\t').build();
                Table targetRegionEnrichment = Table.read().usingOptions(options);
                chromosomePosCodec.addEncodedChrPosColumn(targetRegionEnrichment, true);
                ratioSupplier.setTargetRegionEnrichment(targetRegionEnrichment);
            }

            Table ratios;

            switch (mConfig.mode())
            {
                case TUMOR_ONLY:
                    final Table diploidRegions = new DiploidRegionLoader(mConfig.TumorOnlyDiploidBed, chromosomePosCodec).build();
                    ratios = ratioSupplier.tumorOnly(diploidRegions);
                    break;
                case GERMLIHE_ONLY:
                    ratios = ratioSupplier.germlineOnly();
                    break;
                default:
                    ratios = ratioSupplier.tumorNormalPair();
            }

            final String outputFilename = CobaltRatioFile.generateFilename(
                    mConfig.OutputDir, mConfig.TumorId != null ? mConfig.TumorId : mConfig.ReferenceId);
            CB_LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
            mVersionInfo.write(mConfig.OutputDir);
            CobaltRatioFile.write(outputFilename, ratios.stream().map(r -> rowToCobaltRatio(r, chromosomePosCodec)).collect(Collectors.toList()));

            // do a GC here to free up some memory for ratio segmentation
            System.gc();

            applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);
        }
        finally
        {
            // we must do this to make sure application will exit on exception
            executorService.shutdown();
        }

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        CB_LOGGER.info("Cobalt run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

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

    @NotNull
    public Table loadGCContent(ChromosomePositionCodec chromosomePosCodec) throws IOException
    {
        Table gcProfileTable = Table.create("gcProfiles",
                LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                DoubleColumn.create(CobaltColumns.GC_CONTENT),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        Collection<GCProfile> gcProfileList = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath).values();

        for (GCProfile gcProfile : gcProfileList)
        {
            Row row = gcProfileTable.appendRow();
            long chrPosIndex = chromosomePosCodec.encodeChromosomePosition(gcProfile.chromosome(), gcProfile.start());
            if (chrPosIndex > 0)
            {
                row.setLong(CobaltColumns.ENCODED_CHROMOSOME_POS, chrPosIndex);
            }
            else
            {
                throw new RuntimeException("Unknown chromosome: " + gcProfile.chromosome());
            }
            row.setDouble(CobaltColumns.GC_CONTENT, gcProfile.gcContent());
            row.setBoolean(CobaltColumns.IS_MAPPABLE, gcProfile.isMappable());
            row.setBoolean(CobaltColumns.IS_AUTOSOME, HumanChromosome.fromString(gcProfile.chromosome()).isAutosome());
        }

        return gcProfileTable;
    }
}
