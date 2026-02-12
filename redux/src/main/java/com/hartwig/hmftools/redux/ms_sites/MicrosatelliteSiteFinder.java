package com.hartwig.hmftools.redux.ms_sites;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.loadChrGcProfileMap;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_sites.MsFinderConfig.CHUNK_SIZE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.jitter.MicrosatelliteSite;
import com.hartwig.hmftools.redux.jitter.RefGenomeMicrosatelliteFile;

import org.apache.commons.lang3.Validate;

public class MicrosatelliteSiteFinder
{
    private final MsFinderConfig mConfig;

    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenVersion;
    private final Map<String,List<GCProfile>> mGcProfiles;
    private final Map<String,List<BaseRegion>> mExomeRegions;

    private final Multimap<UnitRepeatKey,MicrosatelliteSite> mAllMsSites;

    public MicrosatelliteSiteFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new MsFinderConfig(configBuilder);
        mGcProfiles = Maps.newHashMap();
        mExomeRegions = Maps.newHashMap();

        mRefGenome = RefGenomeSource.loadRefGenome(mConfig.RefGenomeFile);
        mRefGenVersion = RefGenomeSource.deriveRefGenomeVersion(mRefGenome);

        loadReferenceData();

        mAllMsSites = ArrayListMultimap.create();
    }

    private void loadReferenceData()
    {
        try
        {
            mGcProfiles.putAll(loadChrGcProfileMap(mConfig.GcProfilePath));
        }
        catch(Exception e)
        {
            RD_LOGGER.error("failed to load GC profiles map");
            System.exit(1);
        }

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(mConfig.EnsemblCacheDir, mRefGenVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        for(Map.Entry<String,List<GeneData>> entry : ensemblDataCache.getChrGeneDataMap().entrySet())
        {
            String chromosome = entry.getKey();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            List<BaseRegion> exomeRegions = Lists.newArrayList();
            mExomeRegions.put(chromosome, exomeRegions);

            for(GeneData geneData : entry.getValue())
            {
                TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                if(transcriptData == null)
                    continue;

                for(ExonData exon : transcriptData.exons())
                {
                    exomeRegions.add(new BaseRegion(exon.Start, exon.End));
                }
            }

            BaseRegion.checkMergeOverlaps(exomeRegions, true);
        }
    }

    public void run() throws Exception
    {
        RD_LOGGER.info("finding MS sites from ref genome: {}", mConfig.RefGenomeFile);

        long startTimeMs = System.currentTimeMillis();

        List<ChromosomeSiteFinder> chrSiteFinders = Lists.newArrayList();

        for(Map.Entry<String,Integer> entry : mRefGenome.chromosomeLengths().entrySet())
        {
            String chromosome = entry.getKey();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            chrSiteFinders.add(new ChromosomeSiteFinder(
                    mRefGenome, chromosome, entry.getValue(), CHUNK_SIZE, this::processSiteInfo));
        }

        List<Callable<Void>> threadTasks = chrSiteFinders.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(threadTasks, mConfig.Threads))
        {
            RD_LOGGER.error("error finding sites");
            System.exit(1);
        }

        /*
        MicrosatelliteSiteFinder.findMicrosatellites(mRefGenome, JitterConstants.MIN_MICROSAT_UNIT_COUNT,
            r -> {
                populateMappability(r, mGcProfiles);

                // put all into multimap
                mAllMsSites.put(new UnitRepeatKey(UnitKey.fromUnit(r.unitString()), r.RepeatCount), r);
            });
         */

        SiteDownsampler siteDownsampler = new SiteDownsampler(mConfig, mAllMsSites, mExomeRegions);
        siteDownsampler.downsampleSites();

        // now write only the downsampled sites to file
        String outputFile = RefGenomeMicrosatelliteFile.generateFilename(mConfig.OutputDir, mRefGenVersion);

        Multimap<UnitRepeatKey,MicrosatelliteSite> downSampledmsSites = siteDownsampler.downSampledMsSites();

        try(RefGenomeMicrosatelliteFile refGenomeMicrosatelliteFile = new RefGenomeMicrosatelliteFile(outputFile))
        {
            List<MicrosatelliteSite> sortedMicrosatelliteSites = new ArrayList<>(downSampledmsSites.values());
            sortedMicrosatelliteSites.sort(Comparator.comparing(site -> site.Region));
            sortedMicrosatelliteSites.forEach(refGenomeMicrosatelliteFile::writeRow);
        }

        RD_LOGGER.info("wrote {} microsatellite sites into {}", downSampledmsSites.size(), outputFile);

        RD_LOGGER.info("MS site finding complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processSiteInfo(final MicrosatelliteSite msSite)
    {
        populateMappability(msSite);
        mAllMsSites.put(new UnitRepeatKey(UnitKey.fromUnit(msSite.unitString()), msSite.RepeatCount), msSite);
    }

    private void populateMappability(final MicrosatelliteSite msSite)
    {
        List<GCProfile> gcProfiles = mGcProfiles.get(msSite.chromosome());
        int siteMid = msSite.Region.start() + (msSite.Region.baseLength()) / 2;

        GCProfile endKey = ImmutableGCProfile.builder().from(gcProfiles.get(0)).end(siteMid).build();

        // find the position using binary search
        int index = Collections.binarySearch(gcProfiles, endKey, Comparator.comparingInt(GCProfile::end));

        // we go backwards and forwards
        if(index < 0)
        {
            index = -index - 1;
        }

        if(index < gcProfiles.size())
        {
            GCProfile gcProfile = gcProfiles.get(index);
            Validate.isTrue(gcProfile.overlaps(msSite.Region.genomeRegion()));
            msSite.setMappability(gcProfile.mappablePercentage());
        }
        else
        {
            msSite.setMappability(0);
            RD_LOGGER.warn("microsatellite site({}) gc profile not found", msSite.Region);
        }
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        MsFinderConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MicrosatelliteSiteFinder msSiteFinder = new MicrosatelliteSiteFinder(configBuilder);
        msSiteFinder.run();
    }
}
