package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ListMultimap;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable
{
    // add the AmberConfig options
    @ParametersDelegate
    private final AmberConfig mConfig = new AmberConfig();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    private AmberPersistence mPersistence;
    private VersionInfo mVersionInfo;
    private ImmutableListMultimap<Chromosome,AmberSite> mChromosomeSites;

    public int run() throws IOException, InterruptedException
    {
        mLoggingOptions.setLogLevel();

        mVersionInfo = new VersionInfo("amber.version");
        AMB_LOGGER.info("AMBER version: {}, build timestamp: {}",
                mVersionInfo.version(), mVersionInfo.buildTime().toLocalTime());

        mPersistence = new AmberPersistence(mConfig);

        AMB_LOGGER.info("Loading vcf file {}", mConfig.BafLociPath);
        mChromosomeSites = ImmutableListMultimap.copyOf(AmberSiteFactory.sites(mConfig.BafLociPath));

        if(!mConfig.isValid())
        {
            AMB_LOGGER.error(" invalid config, exiting");
            return 1;
        }

        if(mConfig.isTumorOnly())
        {
            runTumorOnly();
        }
        else if (mConfig.isGermlineOnly())
        {
            runGermlineOnly();
        }
        else
        {
            runNormalMode();
        }
        return 0;
    }

    private void runGermlineOnly() throws InterruptedException, IOException
    {
        AmberGermline germline = new AmberGermline(mConfig, readerFactory(mConfig), mChromosomeSites);

        final List<AmberBAF> amberBAFList = germline.getHeterozygousLoci().values().stream().map(AmberBAF::create)
                .filter(AmberUtils::isValid).sorted().collect(toList());

        mPersistence.persistQC(Collections.emptyList(), germline.getConsanguinityProportion(), germline.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistSnpCheck(germline.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germline.getRegionsOfHomozygosity());
    }

    private void runNormalMode() throws InterruptedException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);

        AmberGermline germline = new AmberGermline(mConfig, readerFactory, mChromosomeSites);

        AmberTumor tumor = new AmberTumor(mConfig, readerFactory,
                germline.getHeterozygousLoci(), germline.getHomozygousLoci());

        final List<TumorBAF> tumorBAFList = tumor.getBafs().values().stream().sorted().collect(toList());
        final List<AmberBAF> amberBAFList = tumorBAFList.stream().map(AmberBAF::create).filter(AmberUtils::isValid).collect(toList());

        final List<TumorContamination> contaminationList = new ArrayList<>(tumor.getContamination().values());

        mPersistence.persistQC(contaminationList, germline.getConsanguinityProportion(), germline.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistContamination(contaminationList);
        mPersistence.persistSnpCheck(germline.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germline.getRegionsOfHomozygosity());
    }

    private void runTumorOnly() throws InterruptedException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);

        final ListMultimap<Chromosome, BaseDepth> allNormal = hetLociTumorOnly();

        // no homozygous sites
        AmberTumor tumor = new AmberTumor(mConfig, readerFactory, allNormal, ArrayListMultimap.create());

        final List<TumorBAF> tumorBAFList = tumor.getBafs().values()
                .stream()
                .filter(x -> x.tumorReadDepth() >= mConfig.TumorOnlyMinDepth)
                .filter(x -> x.tumorRefSupport() >= mConfig.TumorOnlyMinSupport)
                .filter(x -> x.tumorAltSupport() >= mConfig.TumorOnlyMinSupport)
                .filter(x -> isFinite(x.refFrequency()) && Doubles.greaterOrEqual(x.refFrequency(), mConfig.TumorOnlyMinVaf))
                .filter(x -> isFinite(x.altFrequency()) && Doubles.greaterOrEqual(x.altFrequency(), mConfig.TumorOnlyMinVaf))
                .sorted()
                .collect(toList());
        final List<AmberBAF> amberBAFList =
                tumorBAFList.stream().map(AmberBAF::create).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        mPersistence.persistQC(Collections.emptyList(), 0.0, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBAF(amberBAFList);
    }

    // the heterozygous loci snp list that we use contains some regions that could be noisy.
    // this is not a problem if we use those to identify loci that are heterozygous in the
    // germline sample. However, in tumor only mode we would be better off removing those
    // regions.
    private ListMultimap<Chromosome, BaseDepth> hetLociTumorOnly() throws IOException
    {
        List<GenomeRegion> excludedRegions = loadTumorOnlyExcludedSnp();
        final ListMultimap<Chromosome, BaseDepth> result = ArrayListMultimap.create();
        int numBlackListed = 0;

        // filter out everything in loaded genome positions that are in these regions
        for (var entry : mChromosomeSites.entries())
        {
            // check against black list
            boolean blacklisted = false;
            for (GenomeRegion gr : excludedRegions)
            {
                if (gr.contains(entry.getValue()))
                {
                    blacklisted = true;
                    break;
                }
            }
            if (blacklisted)
            {
                numBlackListed++;
            }
            else
            {
                result.put(entry.getKey(), BaseDepthFactory.fromAmberSite(entry.getValue()));
            }
        }
        AMB_LOGGER.info("removed {} blacklisted loci, {} remaining", numBlackListed, result.size());
        return result;
    }

    private static SamReaderFactory readerFactory(final AmberConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.Stringency);
        if(config.RefGenomePath != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomePath)));
        }
        return readerFactory;
    }

    private List<GenomeRegion> loadTumorOnlyExcludedSnp() throws IOException
    {
        String resourcePath = null;
        switch (mConfig.refGenomeVersion)
        {
            case V37:
                // we don't have excluded region for v37 genome
                return Collections.emptyList();
            case V38:
                resourcePath = "tumorOnlyExcludedSnp.38.bed";
                break;
        }

        return AmberUtils.loadBedFromResource(resourcePath);
    }

    @Override
    public void close()
    {
        AMB_LOGGER.info("Complete");
    }

    public static void main(final String... args) throws IOException, InterruptedException
    {
        AmberApplication amberApp = new AmberApplication();
        JCommander commander = JCommander.newBuilder()
                .addObject(amberApp)
                .build();

        // use unix style formatter
        commander.setUsageFormatter(new UnixStyleUsageFormatter(commander));
        // help message show in order parameters are declared
        commander.setParameterDescriptionComparator(new DeclaredOrderParameterComparator(AmberApplication.class));

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
            AMB_LOGGER.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(amberApp.run());
    }

}
