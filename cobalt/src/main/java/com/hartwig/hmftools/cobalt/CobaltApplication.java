package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltColumns.CHROMOSOME;
import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConfig.registerConfig;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.cobalt.CobaltUtils.rowToCobaltRatio;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.calculations.CobaltCalculator;
import com.hartwig.hmftools.cobalt.count.BRC;
import com.hartwig.hmftools.cobalt.count.BamReadCounter;
import com.hartwig.hmftools.cobalt.diploid.DiploidRegionLoader;
import com.hartwig.hmftools.cobalt.exclusions.SuppliedExcludedRegions;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.ColumnType;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.CsvReadOptions;

public class CobaltApplication
{
    private final CobaltConfig mConfig;

    public CobaltApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new CobaltConfig(configBuilder);

        try
        {
            mConfig.validate();
        }
        catch(Exception e)
        {
            CB_LOGGER.error("config loading failed: {}", e.toString());
            System.exit(1);
        }
    }

    private void run()
    {
        long startTimeMs = System.currentTimeMillis();

        CB_LOGGER.info("reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        try
        {
            SamReaderFactory readerFactory = readerFactory(mConfig);

            ChromosomePositionCodec chromosomePosCodec = new ChromosomePositionCodec();

//            BamReadCounter bamReadCounter = new BamReadCounter(WINDOW_SIZE, mConfig, executorService, readerFactory, chromosomePosCodec);
//
//            bamReadCounter.generateDepths();

            BRC brcTumor = null;
            if (mConfig.TumorBamPath != null)
            {
                brcTumor = new BRC(WINDOW_SIZE, mConfig, executorService, mConfig.TumorBamPath, chromosomePosCodec);
            }

            BRC brcReference = null;
            if (mConfig.ReferenceBamPath != null)
            {
                brcReference = new BRC(WINDOW_SIZE, mConfig, executorService, mConfig.ReferenceBamPath, chromosomePosCodec);
            }

            Table tumorReadDepths = brcTumor != null ? brcTumor.generateDepths() : null;
            Table referenceReadDepths = brcReference != null ? brcReference.generateDepths() : null;

            Table gcProfiles = loadMappabilityData(chromosomePosCodec);

            RatioSupplier ratioSupplier = new RatioSupplier(
                    mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir, gcProfiles, referenceReadDepths, tumorReadDepths,
                    chromosomePosCodec);

            if(mConfig.TargetRegionNormFile != null)
            {
                CsvReadOptions options = CsvReadOptions.builder(
                        mConfig.TargetRegionNormFile)
                        .separator(TSV_DELIM.charAt(0))
                        .columnTypesPartial(Map.of(CHROMOSOME, ColumnType.STRING)).build();

                Table targetRegionEnrichment = Table.read().usingOptions(options);

                if(mConfig.SpecificChrRegions.hasFilters())
                {
                    List<ChromosomeData> chromosomeData = brcTumor != null ? brcTumor.chromosomes() : brcReference.chromosomes();
                    List<String> validChromosomes = chromosomeData.stream().map(x -> x.Name).collect(Collectors.toList());

                    targetRegionEnrichment = targetRegionEnrichment.where(
                            targetRegionEnrichment.stringColumn(CobaltColumns.CHROMOSOME).isIn(validChromosomes));
                }

                chromosomePosCodec.addEncodedChrPosColumn(targetRegionEnrichment, true);
                ratioSupplier.setTargetRegionEnrichment(targetRegionEnrichment);
            }

            Table ratios;

            switch(mConfig.mode())
            {
                case TUMOR_ONLY:
                    Table diploidRegions = new DiploidRegionLoader(chromosomePosCodec, mConfig.TumorOnlyDiploidBed).build();
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
            CobaltCalculator calculator = new CobaltCalculator(brcTumor.calculateReadDepths(), mConfig);
            ListMultimap<Chromosome, CobaltRatio> results = calculator.doCalculation();

//            final List<CobaltRatio> collectedRatios = ratios.stream()
//                    .map(r -> rowToCobaltRatio(r, chromosomePosCodec))
//                    .collect(Collectors.toList());
            final List<CobaltRatio> collectedRatios = new ArrayList<>();
            results.keySet().forEach(chromosome -> collectedRatios.addAll(results.get(chromosome)));
            CobaltRatioFile.write(outputFilename, collectedRatios);

            if(!mConfig.SkipPcfCalc)
            {
//                applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);
            }

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

    private static GCProfile mask(GCProfile profile)
    {
        // todo make GProfile a record with an 'excluded' method, or use a new class here
        return ImmutableGCProfile.builder().from(profile).mappablePercentage(0.0).build();
    }

    public Table loadMappabilityData(final ChromosomePositionCodec chromosomePosCodec) throws IOException
    {
        final ListMultimap<Chromosome, GCProfile> gcProfileData = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath);
        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(mConfig.mExcludedRegions);
        ListMultimap<Chromosome, GCProfile> toExclude = excludedRegions.findIntersections(gcProfileData);
        toExclude.forEach((chromosome, gcProfile) ->
                {
                    gcProfileData.remove(chromosome, gcProfile);
                    gcProfileData.put(chromosome, mask(gcProfile));
                }
        );

        Table gcProfileTable = Table.create("gcProfiles",
                LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                DoubleColumn.create("unused"),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));
        Collection<GCProfile> gcProfileList = gcProfileData.values();
        for(GCProfile gcProfile : gcProfileList)
        {
            Row row = gcProfileTable.appendRow();

            String chromosome = gcProfile.chromosome();

            if(mConfig.SpecificChrRegions.hasFilters() && mConfig.SpecificChrRegions.excludeChromosome(chromosome))
            {
                continue;
            }

            long chrPosIndex = chromosomePosCodec.encodeChromosomePosition(chromosome, gcProfile.start());

            if(chrPosIndex > 0)
            {
                row.setLong(CobaltColumns.ENCODED_CHROMOSOME_POS, chrPosIndex);
            }
            else
            {
                throw new RuntimeException("Unknown chromosome: " + chromosome);
            }

            row.setBoolean(CobaltColumns.IS_MAPPABLE, gcProfile.isMappable());
            row.setBoolean(CobaltColumns.IS_AUTOSOME, HumanChromosome.fromString(chromosome).isAutosome());
        }

        return gcProfileTable;
    }

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CobaltApplication application = new CobaltApplication(configBuilder);
        application.run();
    }
}
