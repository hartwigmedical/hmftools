package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.APP_NAME;
import static com.hartwig.hmftools.amber.AmberConstants.TARGET_REGION_SITE_BUFFER;
import static com.hartwig.hmftools.amber.AmberUtils.aboveQualFilter;
import static com.hartwig.hmftools.amber.AmberUtils.fromBaseDepth;
import static com.hartwig.hmftools.amber.AmberUtils.fromTumorBaf;
import static com.hartwig.hmftools.amber.AmberUtils.isValid;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.region.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;

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
import com.hartwig.hmftools.amber.contamination.TumorContamination;
import com.hartwig.hmftools.amber.contamination.TumorContaminationModel;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
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
        ListMultimap<Chromosome, AmberSite> amberSitesMap = AmberSitesFile.sites(mConfig.BafLociPath);

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

                if(amberSites == null)
                {
                    continue;
                }

                int regionIndex = 0;
                BaseRegion currentRegion = regions.get(0);

                for(AmberSite amberSite : amberSites)
                {
                    if(amberSite.position() < currentRegion.start() - TARGET_REGION_SITE_BUFFER)
                    {
                        continue;
                    }

                    while(amberSite.position() > currentRegion.end() + TARGET_REGION_SITE_BUFFER)
                    {
                        ++regionIndex;

                        if(regionIndex >= regions.size())
                        {
                            break;
                        }

                        currentRegion = regions.get(regionIndex);
                    }

                    if(regionIndex >= regions.size())
                    {
                        break;
                    }

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

    private void runGermlineOnly() throws Exception
    {
        GermlineAnalysis germline = new GermlineAnalysis(mConfig, readerFactory(mConfig), mChromosomeSites);

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

        List<TumorBAF> tumorBAFList = tumor.getBafs().values().stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .filter(x -> aboveQualFilter(x.TumorEvidence))
                .sorted().collect(toList());

        List<AmberBAF> amberBAFList = tumorBAFList.stream().map(x -> fromTumorBaf(x)).filter(AmberUtils::isValid).collect(toList());

        List<TumorContamination> contaminationList = new ArrayList<>(tumor.getContamination().values());

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

        ListMultimap<Chromosome, PositionEvidence> allNormal = hetLociTumorOnly();

        // no homozygous sites
        TumorAnalysis tumor = new TumorAnalysis(mConfig, readerFactory, allNormal, ArrayListMultimap.create());

        double contamination = persistRawTumorBAFs(tumor.getBafs());
        List<TumorBAF> tumorBAFList = tumor.getBafs().values()
                .stream()
                .filter(x -> x.TumorEvidence.ReadDepth >= mConfig.TumorMinDepth)
                .filter(x -> aboveQualFilter(x.TumorEvidence))
                .filter(x -> x.TumorEvidence.RefSupport >= mConfig.TumorOnlyMinSupport)
                .filter(x -> x.TumorEvidence.AltSupport >= mConfig.TumorOnlyMinSupport)
                .filter(x -> isFinite(x.refFrequency()) && Doubles.greaterOrEqual(x.refFrequency(), mConfig.TumorOnlyMinVaf))
                .filter(x -> isFinite(x.altFrequency()) && Doubles.greaterOrEqual(x.altFrequency(), mConfig.TumorOnlyMinVaf))
                .sorted()
                .collect(toList());

        List<AmberBAF> amberBAFList = tumorBAFList.stream()
                .map(x -> fromTumorBaf(x)).filter(x -> Double.isFinite(x.tumorBAF())).collect(toList());

        mPersistence.persistQC(0, contamination, null);
        mPersistence.persistVersionInfo(mVersionInfo);
        mPersistence.persistBAF(amberBAFList);
    }

    private double persistRawTumorBAFs(ListMultimap<Chromosome, TumorBAF> tumorBAFs) throws IOException
    {
        List<PositionEvidence> dataList = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(!tumorBAFs.containsKey(chromosome))
            {
                continue;
            }
            List<TumorBAF> chromosomeBAFs = tumorBAFs.get(chromosome);
            dataList.addAll(chromosomeBAFs.stream().map(x -> x.TumorEvidence).toList());
        }

        List<TumorContamination> contaminationList = new ArrayList<>();
        dataList.forEach(x ->
        {
            BaseDepthData bdd = ImmutableBaseDepthData.builder()
                    .ref(switchBases(x.Ref))
                    .alt(switchBases(x.Alt))
                    .refSupport(x.RefSupport)
                    .altSupport(x.AltSupport)
                    .readDepth(x.ReadDepth)
                    .indelCount(x.IndelCount)
                    .build();
            contaminationList.add(new TumorContamination(x.Chromosome, x.Position, null, bdd));
        });

        long sampleHetCount = dataList.size(); // TODO

        double contamination = new TumorContaminationModel().calcContamination(contaminationList, sampleHetCount);

        AMB_LOGGER.info("tumor contamination: {}", contamination);
        mPersistence.persistContamination(contaminationList);

        String filename = PositionEvidenceFile.generateAmberFilenameForWriting(mConfig.OutputDir, mConfig.TumorId);
        PositionEvidenceFile.write(filename, dataList);

        return contamination;
    }

    static BaseDepthData.Base switchBases(PositionEvidence.Base base)
    {
        if(base == PositionEvidence.Base.A)
        {
            return BaseDepthData.Base.A;
        }
        if(base == PositionEvidence.Base.C)
        {
            return BaseDepthData.Base.C;
        }
        if(base == PositionEvidence.Base.G)
        {
            return BaseDepthData.Base.G;
        }
        if(base == PositionEvidence.Base.T)
        {
            return BaseDepthData.Base.T;
        }
        return BaseDepthData.Base.N;
    }

    // the heterozygous loci snp list that we use contains some regions that could be noisy.
    // this is not a problem if we use those to identify loci that are heterozygous in the
    // germline sample. However, in tumor only mode we would be better off removing those regions
    private ListMultimap<Chromosome, PositionEvidence> hetLociTumorOnly() throws IOException
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

        AMB_LOGGER.info("removed {} blacklisted loci, {} remaining", numBlackListed, result.size());
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
