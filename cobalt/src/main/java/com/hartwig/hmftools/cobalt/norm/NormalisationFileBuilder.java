package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.norm.DataLoader.addCobaltSampleData;
import static com.hartwig.hmftools.cobalt.norm.DataLoader.addTargetRegions;
import static com.hartwig.hmftools.cobalt.norm.FileWriter.writeDetailedFile;
import static com.hartwig.hmftools.cobalt.norm.FileWriter.writeNormalisationFile;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MIN_ENRICHMENT_RATIO;
import static com.hartwig.hmftools.cobalt.norm.Normaliser.calcRelativeEnrichment;
import static com.hartwig.hmftools.cobalt.norm.Normaliser.calcSampleAdjustedRatios;
import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.genome.gc.GCBucket.calcGcBucket;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberGender;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class NormalisationFileBuilder
{
    private final NormalisationConfig mConfig;
    private final Map<String,List<RegionData>> mChrRegionData;
    private final GcProfileCache mGcProfileCache;

    public NormalisationFileBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new NormalisationConfig(configBuilder);

        if(mConfig.SampleIds.isEmpty())
        {
            CB_LOGGER.error("no sample IDs loaded");
            System.exit(1);
        }

        mChrRegionData = Maps.newHashMap();

        mGcProfileCache = new GcProfileCache(mConfig.GcProfile);
    }

    public void run()
    {
        CB_LOGGER.info("running Cobalt normalisation file generation from {} samples", mConfig.SampleIds.size());

        // load reference data
        loadTargetRegionsBed(mConfig.TargetRegionsBed);

        setGcProfileData();

        // load Amber files and establish gender
        Map<String,Gender> sampleGenders = !mConfig.SampleGender.isEmpty() ? mConfig.SampleGender : determineAmberGenders();

        // load per-sample Cobalt ratios
        loadSampleCobaltData(sampleGenders);

        // calculate per-sample normalised tumor GC ratios
        calcSampleAdjustedRatios(mConfig.SampleIds, mChrRegionData);

        // calculate a final relative panel enichment ratio for each region
        calcRelativeEnrichment(mChrRegionData, MIN_ENRICHMENT_RATIO);

        writeNormalisationFile(mChrRegionData, mConfig.RefGenVersion, mConfig.OutputFile);

        if(mConfig.DetailedFile != null)
            writeDetailedFile(mChrRegionData, mConfig.SampleIds, mConfig.DetailedFile);

        CB_LOGGER.info("Cobalt normalisation file generation complete");
    }

    private Map<String,Gender> determineAmberGenders()
    {
        Map<String,Gender> sampleGenders = Maps.newHashMap();

        for(String sampleId : mConfig.SampleIds)
        {
            try
            {
                String sampleDir = convertWildcardSamplePath(mConfig.AmberDir, sampleId);
                final String amberFilename = AmberBAFFile.generateAmberFilenameForReading(sampleDir, sampleId);

                Multimap<Chromosome, AmberBAF> chromosomeBafs = AmberBAFFile.read(amberFilename, true);
                Gender gender = AmberGender.determineGender(mConfig.RefGenVersion, chromosomeBafs);
                sampleGenders.put(sampleId, gender);

                CB_LOGGER.debug("sample({}) Amber gender({})", sampleId, gender);
            }
            catch(IOException e)
            {
                CB_LOGGER.error("sample({}) failed to read Amber data: {}", sampleId, e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        return sampleGenders;
    }

    private void loadTargetRegionsBed(final String bedFile)
    {
        if(bedFile == null)
            return;

        List<ChrBaseRegion> regions = ChrBaseRegion.loadChrBaseRegionList(bedFile);

            addTargetRegions(regions, mChrRegionData);

        CB_LOGGER.info("loaded {} target regions from file({})",
                mChrRegionData.values().stream().mapToInt(x -> x.size()).sum(), bedFile);
    }

    private void setGcProfileData()
    {
        for(Map.Entry<String,List<RegionData>> entry : mChrRegionData.entrySet())
        {
            String chromosome = entry.getKey();

            for(RegionData regionData : entry.getValue())
            {
                GCProfile gcProfile = mGcProfileCache.findGcProfile(chromosome, regionData.Position);
                if(gcProfile != null)
                    regionData.setGcProfile(calcGcBucket(gcProfile.gcContent()), gcProfile.mappablePercentage());
            }
        }
    }

    private void loadSampleCobaltData(final Map<String,Gender> sampleGenders)
    {
        for(String sampleId : mConfig.SampleIds)
        {
            Gender amberGender = sampleGenders.get(sampleId);

            String cobaltPanelDir = convertWildcardSamplePath(mConfig.CobaltPanelDir, sampleId);
            String cobaltPanelFilename = CobaltRatioFile.generateFilenameForReading(cobaltPanelDir, sampleId);

            String cobaltWgsFilename = "";

            if(!mConfig.CobaltWgsDir.isEmpty())
            {
                String wgsSampleId = mConfig.getWgsSampleId(sampleId);
                String cobaltWgsDir = convertWildcardSamplePath(mConfig.CobaltWgsDir, wgsSampleId);
                cobaltWgsFilename = CobaltRatioFile.generateFilenameForReading(cobaltWgsDir, wgsSampleId);
            }

            addCobaltSampleData(amberGender, cobaltPanelFilename, cobaltWgsFilename, mChrRegionData);
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        NormalisationConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        NormalisationFileBuilder normFileBuilder = new NormalisationFileBuilder(configBuilder);
        normFileBuilder.run();
    }
}
