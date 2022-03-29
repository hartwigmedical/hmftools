package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.cobalt.RatioSegmentation.applyRatioSegmentation;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class CobaltApplication implements AutoCloseable
{
    @ParametersDelegate
    private final CobaltConfig mConfig = new CobaltConfig();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
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

        if (mConfig.RefGenomePath != null && !(new File(mConfig.RefGenomePath).exists()))
        {
            throw new IOException("Unable to locate ref genome file " + mConfig.RefGenomePath);
        }

        CB_LOGGER.info("Reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        final Multimap<Chromosome, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath);

        final List<BEDFeature> diploidBedFile = diploidBedFile(mConfig);

        final SamReaderFactory readerFactory = readerFactory(mConfig);

        final CountSupplier countSupplier = new CountSupplier(
                WINDOW_SIZE, mConfig.MinMappingQuality,
                executorService, readerFactory);
        
        final Multimap<Chromosome, CobaltCount> readCounts = countSupplier.generateCounts(mConfig.ReferenceBamPath, mConfig.TumorBamPath);

        final RatioSupplier ratioSupplier = new RatioSupplier(mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir);

        Multimap<Chromosome, CobaltRatio> ratios;

        switch (mConfig.mode())
        {
            case TUMOR_ONLY:
                ratios = ratioSupplier.tumorOnly(diploidBedFile, gcProfiles, readCounts);
                break;
            case GERMLIHE_ONLY:
                ratios = ratioSupplier.germlineOnly(gcProfiles, readCounts);
                break;
            default:
                ratios = ratioSupplier.tumorNormalPair(gcProfiles, readCounts);
        }

        final String outputFilename = CobaltRatioFile.generateFilenameForWriting(
                mConfig.OutputDir, mConfig.TumorId != null ? mConfig.TumorId : mConfig.ReferenceId);
        CB_LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
        mVersionInfo.write(mConfig.OutputDir);
        CobaltRatioFile.write(outputFilename, ratios);

        applyRatioSegmentation(executorService, mConfig.OutputDir, outputFilename, mConfig.ReferenceId, mConfig.TumorId);

        executorService.shutdown();

        CB_LOGGER.info("Complete");

        return 0;
    }

    @NotNull
    private static List<BEDFeature> diploidBedFile(CobaltConfig config) throws IOException
    {
        List<BEDFeature> result = new ArrayList<>();
        if(config.mode() != CobaltConfig.Mode.TUMOR_ONLY)
        {
            return result;
        }

        CB_LOGGER.info("Reading diploid regions from {}", config.TumorOnlyDiploidBed);
        try(final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(config.TumorOnlyDiploidBed,
                new BEDCodec(),
                false))
        {
            for(BEDFeature bedFeature : reader.iterator())
            {
                result.add(bedFeature);
            }
        }

        return result;
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
}
