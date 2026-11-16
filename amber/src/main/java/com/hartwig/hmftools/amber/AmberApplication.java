package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.lang.String.format;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.APP_NAME;
import static com.hartwig.hmftools.amber.AmberConstants.TARGET_REGION_SITE_BUFFER;
import static com.hartwig.hmftools.amber.AmberUtils.aboveNoiseFloor;
import static com.hartwig.hmftools.amber.AmberUtils.aboveQualFilter;
import static com.hartwig.hmftools.amber.AmberUtils.aboveTumorOnlyMinSupport;
import static com.hartwig.hmftools.amber.AmberUtils.fromBaseDepth;
import static com.hartwig.hmftools.amber.AmberUtils.isValid;
import static com.hartwig.hmftools.amber.LohIdentifier.checkGenderSites;
import static com.hartwig.hmftools.common.amber.AmberGender.determineGender;
import static com.hartwig.hmftools.common.region.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ListMultimap;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.amber.contamination.TumorContamination;
import com.hartwig.hmftools.amber.contamination.TumorContaminationModel;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.amber.blacklist.AmberBlacklistFile;
import com.hartwig.hmftools.amber.blacklist.AmberBlacklistPoint;
import com.hartwig.hmftools.amber.purity.PurityAnalysisConfig;
import com.hartwig.hmftools.amber.purity.TumorOnlyPurityAnalysis;
import com.hartwig.hmftools.amber.purity.CandidatePeak;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.VersionInfo;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class AmberApplication implements AutoCloseable
{
    private final AmberConfig mConfig;

    private ResultsWriter mPersistence;
    private VersionInfo mVersionInfo;
    private ImmutableListMultimap<Chromosome, AmberSite> mChromosomeSites;

    public AmberApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new AmberConfig(configBuilder);
    }

    public int run() throws Exception
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

    private ImmutableListMultimap<Chromosome, AmberSite> loadAmberSites() throws IOException
    {
        ListMultimap<Chromosome,AmberSite> amberSitesMap = AmberSitesFile.sites(mConfig.BafLociPath);

        if(mConfig.SpecificChrRegions.hasFilters())
        {
            ListMultimap<Chromosome,AmberSite> filteredSitesMap = ArrayListMultimap.create();

            for(AmberSite amberSite : amberSitesMap.values())
            {
                if(mConfig.SpecificChrRegions.includePosition(amberSite.Chromosome, amberSite.Position))
                    filteredSitesMap.put(amberSite.chr(), amberSite);
            }

            amberSitesMap = filteredSitesMap;
        }

        if(mConfig.TargetRegionsBed == null)
        {
            return ImmutableListMultimap.copyOf(amberSitesMap);
        }

        ListMultimap<Chromosome, AmberSite> targetRegionSites = ArrayListMultimap.create();
        try
        {
            Map<Chromosome, List<BaseRegion>> targetRegions = loadBedFileChrMap(mConfig.TargetRegionsBed);

            for(Map.Entry<Chromosome, List<BaseRegion>> entry : targetRegions.entrySet())
            {
                Chromosome chromosome = entry.getKey();
                List<BaseRegion> regions = entry.getValue();

                Collection<AmberSite> amberSites = amberSitesMap.get(chromosome);

                int regionIndex = 0;
                BaseRegion currentRegion = regions.get(0);

                for(AmberSite amberSite : amberSites)
                {
                    if(amberSite.position() < currentRegion.start() - TARGET_REGION_SITE_BUFFER)
                        continue;

                    while(amberSite.position() > currentRegion.end() + TARGET_REGION_SITE_BUFFER)
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
            AMB_LOGGER.error(format("failed to load target regions file(): {%s}", mConfig.TargetRegionsBed), e);
            System.exit(1);
        }

        return ImmutableListMultimap.copyOf(targetRegionSites);
    }

    private void runGermlineOnly() throws Exception
    {
        GermlineAnalysis germline = new GermlineAnalysis(mConfig, readerFactory(mConfig), mChromosomeSites);

        AMB_LOGGER.debug("BAM reading complete");

        List<AmberBAF> amberBAFList = Lists.newArrayList();

        for(PositionEvidence baseDepth : germline.getHeterozygousLoci().values())
        {
            AmberBAF amberBAF = fromBaseDepth(baseDepth);

            if(mConfig.WriteUnfilteredGermline || isValid(amberBAF))
            {
                amberBAFList.add(amberBAF);
            }
        }

        Collections.sort(amberBAFList);

        mPersistence.persistQC(germline.getConsanguinityProportion(), 0, germline.getUniparentalDisomy());

        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistSnpCheck(germline.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germline.getRegionsOfHomozygosity());
    }

    private void runNormalMode() throws Exception
    {
        SamReaderFactory readerFactory = readerFactory(mConfig);

        GermlineAnalysis germline = new GermlineAnalysis(mConfig, readerFactory, mChromosomeSites);

        TumorAnalysis tumor = new TumorAnalysis(mConfig, readerFactory, germline.getHeterozygousLoci(), germline.getHomozygousLoci());

        AMB_LOGGER.debug("BAM reading complete");

        List<TumorBAF> tumorBAFList = tumor.chrBafMap().values().stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .filter(x -> aboveQualFilter(x.TumorEvidence))
                .sorted().toList();

        List<AmberBAF> amberBAFList = tumorBAFList.stream().map(AmberUtils::fromTumorBaf).filter(AmberUtils::isValid).collect(toList());

        if(mConfig.WriteTumorData)
        {
            List<PositionEvidence> rawData = tumorBAFList.stream().map(x -> x.TumorEvidence).toList();
            String rawDataFileName = PositionEvidenceFile.generateTumorDataFilename(mConfig.OutputDir, mConfig.TumorId);
            PositionEvidenceFile.write(rawDataFileName, rawData);
        }

        List<TumorContamination> contaminationList = tumor.contaminationSites();

        long sampleHetCount = amberBAFList.size();

        double contamination = new TumorContaminationModel().calcContamination(contaminationList, sampleHetCount);

        mPersistence.persistQC(germline.getConsanguinityProportion(), contamination, germline.getUniparentalDisomy());
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistContamination(contaminationList);
        mPersistence.persistSnpCheck(germline.getSnpCheckedLoci());
        mPersistence.persistBAF(amberBAFList);
        mPersistence.persistHomozygousRegions(germline.getRegionsOfHomozygosity());
    }

    private void runTumorOnly() throws Exception
    {
        SamReaderFactory readerFactory = readerFactory(mConfig);

        ListMultimap<Chromosome,PositionEvidence> allNormal = hetLociTumorOnly();

        // no homozygous sites
        TumorAnalysis tumor = new TumorAnalysis(mConfig, readerFactory, allNormal, ArrayListMultimap.create());

        AMB_LOGGER.debug("BAM reading complete");

        LohIdentifier lohIdentifier = new LohIdentifier(mConfig.RefGenVersion, mConfig.TumorMinDepth);

        lohIdentifier.findLohSegments(tumor.chrBafMap());
        AMB_LOGGER.debug("LOH identification complete");

        List<TumorBAF> readDepthAndQualityFiltered = tumor.chrBafMap().values()
                .stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .filter(x -> aboveQualFilter(x.TumorEvidence))
                .filter(x -> aboveTumorOnlyMinSupport(x.TumorEvidence, mConfig.TumorOnlyMinSupport) || x.TumorEvidence.lohRecovered())
                .sorted().toList();

        List<PositionEvidence> rawData = readDepthAndQualityFiltered.stream().map(x -> x.TumorEvidence).toList();

        if(mConfig.WriteTumorData)
        {
            String rawDataFileName = PositionEvidenceFile.generateTumorDataFilename(mConfig.OutputDir, mConfig.TumorId);
            PositionEvidenceFile.write(rawDataFileName, rawData);
        }

        TumorOnlyPurityAnalysis noiseFloorAnalysis = new TumorOnlyPurityAnalysis(
                rawData.stream().filter(x -> !x.lohRecovered()).collect(toList()), // ignore recovered points
                mChromosomeSites, new PurityAnalysisConfig(mConfig));

        double noiseFloor = noiseFloorAnalysis.cutoff();

        double contamination = noiseFloorAnalysis.contaminationPeaks().stream().map(CandidatePeak::vaf).max(Double::compare).orElse(0.0);
        AMB_LOGGER.debug(format("contamination level(%.3f) noiseFloor(%.3f)", contamination, noiseFloor));

        List<TumorBAF> tumorBAFList = readDepthAndQualityFiltered
                .stream()
                .filter(x -> aboveNoiseFloor(x, noiseFloor))
                .sorted().toList();

        List<AmberBAF> amberBAFList = tumorBAFList.stream()
                .map(AmberUtils::fromTumorBaf).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        // amberBAFList = checkGenderSites(mConfig.RefGenVersion, amberBAFList); // no need to check chr-X is not handled by LOH identifier

        mPersistence.persistQC(0, contamination, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBAF(amberBAFList);
    }

    // the heterozygous loci snp list that we use contains some regions that could be noisy.
    // this is not a problem if we use those to identify loci that are heterozygous in the
    // germline sample. However, in tumor-only mode we would be better off removing those regions
    private ListMultimap<Chromosome, PositionEvidence> hetLociTumorOnly()
    {
        List<GenomeRegion> excludedRegions = loadTumorOnlyExcludedSnp();
        ListMultimap<Chromosome, PositionEvidence> result = ArrayListMultimap.create();
        int numBlackListed = 0;

        // filter out everything in loaded genome positions that are in these regions
        for(Map.Entry<Chromosome, AmberSite> entry : mChromosomeSites.entries())
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

        if(numBlackListed > 0)
        {
            AMB_LOGGER.info("removed {} blacklisted loci, {} remaining", numBlackListed, result.size());
        }

        return result;
    }

    private static SamReaderFactory readerFactory(final AmberConfig config)
    {
        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }
        return readerFactory;
    }

    private List<GenomeRegion> loadTumorOnlyExcludedSnp()
    {
        if(mConfig.RefGenVersion == V37)
        {
            return Collections.emptyList();
        }

        List<ChrBaseRegion> regions = AmberUtils.loadBedFromResource("tumorOnlyExcludedSnp.38.bed");
        return regions.stream().map(ChrBaseRegion::genomeRegion).collect(toList());
    }

    @Override
    public void close()
    {
        AMB_LOGGER.info("Amber complete");
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AmberConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AmberApplication amberApp = new AmberApplication(configBuilder);
        amberApp.run();
    }
}
