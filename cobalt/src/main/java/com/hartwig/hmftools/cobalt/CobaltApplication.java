package com.hartwig.hmftools.cobalt;

import static com.google.common.collect.ArrayListMultimap.create;
import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConfig.registerConfig;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;

import java.io.IOException;
import java.util.ArrayList;
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
import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.CobaltMedianRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;

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

            if(mConfig.TumorId != null)
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
                applyRatioSegmentation(executorService, mConfig.OutputDir, mConfig.cobaltRatiosFileName(), mConfig.ReferenceId, mConfig.TumorId, mConfig.PcfGamma);
            }

            final VersionInfo version = fromAppName(APP_NAME);
            version.write(mConfig.OutputDir);
        }
        catch(Exception e)
        {
            CB_LOGGER.error("error running Cobalt", e);
            System.exit(1);
        }
        finally
        {
            executorService.shutdown();
        }

        CB_LOGGER.info("Cobalt complete, mins({})", runTimeMinsStr(startTimeMs));
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
