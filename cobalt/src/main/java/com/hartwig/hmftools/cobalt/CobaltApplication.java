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
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class CobaltApplication implements AutoCloseable
{
    private final CobaltConfig mConfig;
    private final VersionInfo mVersionInfo;

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException
    {
        final Options options = CobaltConfig.createOptions();
        try(final CobaltApplication application = new CobaltApplication(options, args))
        {
            application.run();
        } 
        catch(ParseException e)
        {
            CB_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CobaltApplication", options);
            System.exit(1);
        }
    }

    private CobaltApplication(final Options options, final String... args) throws ParseException, IOException
    {
        mVersionInfo = new VersionInfo("cobalt.version");
        CB_LOGGER.info("COBALT version: {}", mVersionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        mConfig = new CobaltConfig(cmd);

        if(mConfig.ReferenceBamPath != null && !mConfig.RefGenomePath.isEmpty() && !new File(mConfig.GcProfilePath).exists())
        {
            throw new IOException("Unable to locate ref genome file " + mConfig.RefGenomePath);
        }
    }

    private void run() throws IOException, ExecutionException, InterruptedException
    {
        if(!mConfig.isValid())
        {
            CB_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        CB_LOGGER.info("Reading GC Profile from {}", mConfig.GcProfilePath);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.ThreadCount, namedThreadFactory);

        final Multimap<Chromosome, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(WINDOW_SIZE, mConfig.GcProfilePath);

        final List<BEDFeature> diploidBedFile = diploidBedFile(mConfig);

        final SamReaderFactory readerFactory = readerFactory(mConfig);

        final CountSupplier countSupplier = new CountSupplier(
                mConfig.TumorId, mConfig.OutputDir, WINDOW_SIZE, mConfig.MinMappingQuality,
                executorService, readerFactory);
        
        final Multimap<Chromosome, CobaltCount> readCounts = mConfig.TumorOnly
                ? countSupplier.tumorOnly(mConfig.TumorBamPath)
                : countSupplier.pairedTumorNormal(mConfig.ReferenceBamPath, mConfig.TumorBamPath);

        final RatioSupplier ratioSupplier = new RatioSupplier(mConfig.ReferenceId, mConfig.TumorId, mConfig.OutputDir);

        final Multimap<Chromosome, CobaltRatio> ratios = mConfig.TumorOnly
                ? ratioSupplier.tumorOnly(diploidBedFile, gcProfiles, readCounts)
                : ratioSupplier.tumorNormalPair(gcProfiles, readCounts);

        final String outputFilename = CobaltRatioFile.generateFilenameForWriting(mConfig.OutputDir, mConfig.TumorId);
        CB_LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
        mVersionInfo.write(mConfig.OutputDir);
        CobaltRatioFile.write(outputFilename, ratios);

        applyRatioSegmentation(executorService, mConfig.OutputDir, mConfig.ReferenceId, mConfig.TumorId);

        executorService.shutdown();
    }

    @NotNull
    private static List<BEDFeature> diploidBedFile(CobaltConfig config) throws IOException
    {
        List<BEDFeature> result = new ArrayList<>();
        if(!config.TumorOnly)
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

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close()
    {
        CB_LOGGER.info("Complete");
    }
}
