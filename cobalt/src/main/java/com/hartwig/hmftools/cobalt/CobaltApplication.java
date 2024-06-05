package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConfig.registerConfig;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.cobalt.CobaltUtils.rowToCobaltRatio;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.count.BamReadCounter;
import com.hartwig.hmftools.cobalt.diploid.DiploidRegionLoader;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import tech.tablesaw.api.*;
import tech.tablesaw.io.csv.CsvReadOptions;

public class CobaltApplication
{
    private final CobaltConfig mConfig;

    public CobaltApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new CobaltConfig(configBuilder);

        try{
            mConfig.validate();
        }
        catch(Exception e)
        {
            CB_LOGGER.error("config loading failed: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CobaltApplication application = new CobaltApplication(configBuilder);
        application.run();
    }

    private void run()
    {
        long startTimeMs = System.currentTimeMillis();

        CB_LOGGER.info("reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        try
        {
            final SamReaderFactory readerFactory = readerFactory(mConfig);

            ChromosomePositionCodec chromosomePosCodec = new ChromosomePositionCodec();

            final BamReadCounter bamReadCounter = new BamReadCounter(WINDOW_SIZE, mConfig, executorService, readerFactory, chromosomePosCodec);

            bamReadCounter.generateDepths(mConfig.ReferenceBamPath, mConfig.TumorBamPath);

            Table referenceReadDepths = bamReadCounter.getReferenceDepths();
            Table tumorReadDepths = bamReadCounter.getTumorDepths();

            final Table gcProfiles = loadGCContent(chromosomePosCodec);

            final RatioSupplier ratioSupplier = new RatioSupplier(mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir,
                    gcProfiles, referenceReadDepths, tumorReadDepths,
                    chromosomePosCodec);

            if(mConfig.TargetRegionPath != null)
            {
                CsvReadOptions options = CsvReadOptions.builder(mConfig.TargetRegionPath)
                        .separator(TSV_DELIM.charAt(0))
                        .columnTypesPartial(Map.of("chromosome", ColumnType.STRING)).build();
                Table targetRegionEnrichment = Table.read().usingOptions(options);
                chromosomePosCodec.addEncodedChrPosColumn(targetRegionEnrichment, true);
                ratioSupplier.setTargetRegionEnrichment(targetRegionEnrichment);
            }

            Table ratios;

            switch(mConfig.mode())
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

            CB_LOGGER.info("persisting cobalt ratios to {}", outputFilename);


            CobaltRatioFile.write(outputFilename, ratios.stream().map(r -> rowToCobaltRatio(r, chromosomePosCodec)).collect(Collectors.toList()));

            applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);

            final VersionInfo version = fromAppName(APP_NAME);
            version.write(mConfig.OutputDir);
        }
        catch(Exception e)
        {
            CB_LOGGER.error("error: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        finally
        {
            executorService.shutdown();
        }

        CB_LOGGER.info("Cobalt complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static SamReaderFactory readerFactory(@NotNull final CobaltConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);

        if(config.RefGenomePath != null && !config.RefGenomePath.isEmpty())
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomePath)));
        }

        return readerFactory;
    }

    public Table loadGCContent(ChromosomePositionCodec chromosomePosCodec) throws IOException
    {
        Table gcProfileTable = Table.create("gcProfiles",
                LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                DoubleColumn.create(CobaltColumns.GC_CONTENT),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        Collection<GCProfile> gcProfileList = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath).values();

        for(GCProfile gcProfile : gcProfileList)
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
