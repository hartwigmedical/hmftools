package com.hartwig.hmftools.cobalt;

import static com.google.common.collect.ArrayListMultimap.create;
import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConfig.registerConfig;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.ListMultimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.calculations.CobaltCalculator;
import com.hartwig.hmftools.cobalt.count.BRC;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.exclusions.SuppliedExcludedRegions;
import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.CobaltMedianRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
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
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;

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
            BRC brcTumor = mConfig.tumorBamReader(executorService);
            BRC brcRef = mConfig.referenceBamReader(executorService);

            ListMultimap<HumanChromosome, DepthReading> tumourDepths = brcTumor == null ? create() : brcTumor.calculateReadDepths();
            CB_LOGGER.info("Tumor depths collected, size: {}", tumourDepths.size());
            ListMultimap<HumanChromosome, DepthReading> refDepths = brcRef == null ? create() : brcRef.calculateReadDepths();
            CB_LOGGER.info("Reference depths collected, size: {}", refDepths.size());

            CobaltCalculator calculator = new CobaltCalculator(tumourDepths, refDepths, mConfig);
            ListMultimap<Chromosome, CobaltRatio> results = calculator.getCalculatedRatios();

            final List<CobaltRatio> collectedRatios = new ArrayList<>();
            results.keySet().forEach(chromosome -> collectedRatios.addAll(results.get(chromosome)));
            CobaltRatioFile.write(mConfig.cobaltRatiosFileName(), collectedRatios);

            if (mConfig.TumorId != null)
            {
                CB_LOGGER.info("persisting tumor {} GC read count to {}", mConfig.TumorId, mConfig.OutputDir);
                CobaltGcMedianFile.write(mConfig.tumorGcMedianFileName(), calculator.tumorMedianReadDepth());
            }

            if(mConfig.ReferenceId != null)
            {
                CB_LOGGER.info("persisting {} gc ratio medians to {}", mConfig.ReferenceId, mConfig.OutputDir);
                CobaltMedianRatioFile.write(mConfig.medianRatiosFileName(), calculator.medianRatios());

                CB_LOGGER.info("persisting reference {} GC read count to {}", mConfig.ReferenceId, mConfig.OutputDir);
                CobaltGcMedianFile.write(mConfig.referenceGcMedianFileName(), calculator.referenceMedianReadDepth());
            }

            if(!mConfig.SkipPcfCalc)
            {
                //                applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);
            }

            final VersionInfo version = fromAppName(APP_NAME);
            version.write(mConfig.OutputDir);
            //     List<DepthReading> tumorRawReads = CobaltUtils.toList(tumourDepths);
            //            final String tumorRawReadsFile = DepthReadingsFile.generateFilename(mConfig.OutputDir, mConfig.TumorId, true);
            //            DepthReadingsFile.write(tumorRawReadsFile, tumorRawReads);
            //            CB_LOGGER.info("Tumor raw reads persisted to: {}", tumorRawReadsFile);
            //
            //            List<DepthReading> refRawReads = CobaltUtils.toList(refDepths);
            //            final String refRawReadsFile = DepthReadingsFile.generateFilename(mConfig.OutputDir, mConfig.ReferenceId, false);
            //            DepthReadingsFile.write(refRawReadsFile, refRawReads);
            //            CB_LOGGER.info("Reference raw reads persisted to: {}", refRawReadsFile);
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
