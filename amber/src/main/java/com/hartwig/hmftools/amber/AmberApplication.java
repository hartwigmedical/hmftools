package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.APP_NAME;
import static com.hartwig.hmftools.amber.AmberConstants.TARGET_REGION_SITE_BUFFER;
import static com.hartwig.hmftools.amber.AmberUtils.fromBaseDepth;
import static com.hartwig.hmftools.amber.AmberUtils.fromTumorBaf;
import static com.hartwig.hmftools.amber.AmberUtils.isValid;
import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ListMultimap;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable
{
    private final AmberConfig mConfig;

    private ResultsWriter mPersistence;
    private VersionInfo mVersionInfo;
    private ImmutableListMultimap<Chromosome,AmberSite> mChromosomeSites;

    public AmberApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new AmberConfig(configBuilder);
    }

    public int run() throws IOException, InterruptedException
    {
        long startTimeMs = System.currentTimeMillis();

        mVersionInfo = fromAppName(APP_NAME);

        mPersistence = new ResultsWriter(mConfig);

        mChromosomeSites = loadAmberSites();

        if(!mConfig.isValid())
        {
            AMB_LOGGER.error(" invalid config, exiting");
            return 1;
        }

        if(mConfig.isTumorOnly())
        {
            runTumorOnly();
        }
        else if(mConfig.isGermlineOnly())
        {
            runGermlineOnly();
        }
        else
        {
            runNormalMode();
        }

        AMB_LOGGER.info("Amber complete, mins({})", runTimeMinsStr(startTimeMs));

        return 0;
    }

    private ImmutableListMultimap<Chromosome,AmberSite> loadAmberSites() throws IOException
    {
        ListMultimap<Chromosome,AmberSite> amberSitesMap = AmberSitesFile.sites(mConfig.BafLociPath);

        if(mConfig.TargetRegionsBed == null)
            return ImmutableListMultimap.copyOf(amberSitesMap);

        ListMultimap<Chromosome,AmberSite> targetRegionSites = ArrayListMultimap.create();

        try
        {
            Map<Chromosome,List<BaseRegion>> targetRegions = loadBedFileChrMap(mConfig.TargetRegionsBed);

            for(Map.Entry<Chromosome,List<BaseRegion>> entry : targetRegions.entrySet())
            {
                Chromosome chromosome = entry.getKey();
                List<BaseRegion> regions = entry.getValue();

                Collection<AmberSite> amberSites = amberSitesMap.get(chromosome);

                if(amberSites == null)
                    continue;

                int regionIndex = 0;
                BaseRegion currentRegion = regions.get(0);

                for(AmberSite amberSite : amberSites)
                {
                    if(amberSite.position() < currentRegion.start())
                        continue;

                    while(amberSite.position() > currentRegion.end())
                    {
                        ++regionIndex;

                        if(regionIndex >= regions.size())
                            break;

                        currentRegion = regions.get(regionIndex);
                    }

                    if(regionIndex >= regions.size())
                        break;

                    if(amberSite.position() >= currentRegion.start() - TARGET_REGION_SITE_BUFFER
                    && amberSite.position() <= currentRegion.end() + TARGET_REGION_SITE_BUFFER)
                    {
                        targetRegionSites.put(chromosome, amberSite);
                    }
                }
            }

            return ImmutableListMultimap.copyOf(targetRegionSites);
        }
        catch(Exception e)
        {
            AMB_LOGGER.error("failed to load target regions file(): {}", mConfig.TargetRegionsBed, e.toString());
            System.exit(1);
        }

        return ImmutableListMultimap.copyOf(targetRegionSites);
    }

    private void runGermlineOnly() throws InterruptedException, IOException
    {
        GermlineAnalysis germline = new GermlineAnalysis(mConfig, readerFactory(mConfig), mChromosomeSites);

        List<AmberBAF> amberBAFList = Lists.newArrayList();

        for(PositionEvidence baseDepth : germline.getHeterozygousLoci().values())
        {
            AmberBAF amberBAF = fromBaseDepth(baseDepth);

            if(mConfig.WriteUnfilteredGermline || isValid(amberBAF))
                amberBAFList.add(amberBAF);
        }

        Collections.sort(amberBAFList);

        mPersistence.persistQC(Collections.emptyList(), germline.getConsanguinityProportion(), germline.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistSnpCheck(germline.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germline.getRegionsOfHomozygosity());
    }

    private void runNormalMode() throws InterruptedException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mConfig);

        GermlineAnalysis germline = new GermlineAnalysis(mConfig, readerFactory, mChromosomeSites);

        TumorAnalysis tumor = new TumorAnalysis(mConfig, readerFactory,
                germline.getHeterozygousLoci(), germline.getHomozygousLoci());

        final List<TumorBAF> tumorBAFList = tumor.getBafs().values().stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .sorted().collect(toList());

        final List<AmberBAF> amberBAFList = tumorBAFList.stream().map(x -> fromTumorBaf(x)).filter(AmberUtils::isValid).collect(toList());

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

        final ListMultimap<Chromosome, PositionEvidence> allNormal = hetLociTumorOnly();

        // no homozygous sites
        TumorAnalysis tumor = new TumorAnalysis(mConfig, readerFactory, allNormal, ArrayListMultimap.create());

        final List<TumorBAF> tumorBAFList = tumor.getBafs().values()
                .stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .filter(x -> x.TumorEvidence.RefSupport >= mConfig.TumorOnlyMinSupport)
                .filter(x -> x.TumorEvidence.AltSupport >= mConfig.TumorOnlyMinSupport)
                .filter(x -> isFinite(x.refFrequency()) && Doubles.greaterOrEqual(x.refFrequency(), mConfig.TumorOnlyMinVaf))
                .filter(x -> isFinite(x.altFrequency()) && Doubles.greaterOrEqual(x.altFrequency(), mConfig.TumorOnlyMinVaf))
                .sorted()
                .collect(toList());

        final List<AmberBAF> amberBAFList = tumorBAFList.stream()
                .map(x -> fromTumorBaf(x)).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        mPersistence.persistQC(Collections.emptyList(), 0.0, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBAF(amberBAFList);
    }

    // the heterozygous loci snp list that we use contains some regions that could be noisy.
    // this is not a problem if we use those to identify loci that are heterozygous in the
    // germline sample. However, in tumor only mode we would be better off removing those
    // regions.
    private ListMultimap<Chromosome,PositionEvidence> hetLociTumorOnly() throws IOException
    {
        List<GenomeRegion> excludedRegions = loadTumorOnlyExcludedSnp();
        final ListMultimap<Chromosome, PositionEvidence> result = ArrayListMultimap.create();
        int numBlackListed = 0;

        // filter out everything in loaded genome positions that are in these regions
        for(Map.Entry<Chromosome,AmberSite> entry : mChromosomeSites.entries())
        {
            // check against black list
            boolean blacklisted = false;
            for(GenomeRegion gr : excludedRegions)
            {
                if(gr.contains(entry.getValue()))
                {
                    blacklisted = true;
                    break;
                }
            }
            if(blacklisted)
            {
                numBlackListed++;
            }
            else
            {
                result.put(entry.getKey(), PositionEvidenceChecker.fromAmberSite(entry.getValue()));
            }
        }

        AMB_LOGGER.info("removed {} blacklisted loci, {} remaining", numBlackListed, result.size());
        return result;
    }

    private static SamReaderFactory readerFactory(final AmberConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }
        return readerFactory;
    }

    private List<GenomeRegion> loadTumorOnlyExcludedSnp()
    {
        if(mConfig.RefGenVersion == V37)
            return Collections.emptyList();

        List<ChrBaseRegion> regions = AmberUtils.loadBedFromResource("tumorOnlyExcludedSnp.38.bed");
        return regions.stream().map(x -> x.genomeRegion()).collect(toList());
    }

    @Override
    public void close()
    {
        AMB_LOGGER.info("Amber complete");
    }

    public static void main(final String... args) throws IOException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AmberConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AmberApplication amberApp = new AmberApplication(configBuilder);
        amberApp.run();
    }
}
