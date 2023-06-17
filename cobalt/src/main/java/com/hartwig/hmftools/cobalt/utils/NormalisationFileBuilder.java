package com.hartwig.hmftools.cobalt.utils;

import static java.lang.Math.floor;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.utils.CobaltDataLoader.addCobaltSampleData;
import static com.hartwig.hmftools.cobalt.utils.NormConstants.REGION_SIZE;
import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.genome.gc.GCBucket.calcGcBucket;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberGender;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class NormalisationFileBuilder
{
    private final NormalisationConfig mConfig;
    private final Map<String,List<RegionData>> mChrRegionData;
    private final GcProfileCache mGcProfileCache;

    public NormalisationFileBuilder(final CommandLine cmd)
    {
        mConfig = new NormalisationConfig(cmd);

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

        Map<String,Gender> sampleGenders = determineAmberGenders();

        loadSampleCobaltData(sampleGenders);

        CB_LOGGER.info("Cobalt normalisation file generation complete");
    }

    private Map<String,Gender> determineAmberGenders()
    {
        Map<String,Gender> sampleGenders = Maps.newHashMap();

        for(String sampleId : mConfig.SampleIds)
        {
            try
            {
                String sampleDir = mConfig.AmberDir.replaceAll("\\*", sampleId);
                final String amberFilename = AmberBAFFile.generateAmberFilenameForReading(sampleDir, sampleId);

                Multimap<Chromosome, AmberBAF> chromosomeBafs = AmberBAFFile.read(amberFilename, true);
                Gender gender = AmberGender.determineGender(mConfig.RefGenVersion, chromosomeBafs);
                sampleGenders.put(sampleId, gender);

                CB_LOGGER.debug("AMBER_GENDER:{},{}", sampleId, gender);
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

        try
        {
            List<NamedBed> namedBedRecords = readBedFile(bedFile);

            String currentChromosome = "";
            List<RegionData> regions = null;

            for(NamedBed namedBed : namedBedRecords)
            {
                if(!namedBed.chromosome().equals(currentChromosome))
                {
                    currentChromosome = namedBed.chromosome();
                    regions = Lists.newArrayList();
                    mChrRegionData.put(namedBed.chromosome(), regions);
                }

                int startPosition = (int)(floor(namedBed.start()/(double)REGION_SIZE) * REGION_SIZE + 1);
                int endPosition = startPosition + REGION_SIZE - 1;

                regions.add(new RegionData(startPosition));

                while(endPosition < namedBed.end())
                {
                    startPosition += REGION_SIZE;
                    regions.add(new RegionData(startPosition));
                    endPosition = startPosition + REGION_SIZE - 1;
                }
            }

            CB_LOGGER.info("loaded {} target regions from file({})",
                    mChrRegionData.values().stream().mapToInt(x -> x.size()).sum(), bedFile);
        }
        catch (IOException e)
        {
            CB_LOGGER.error("failed to load target regions BED file: {}", e.toString());
        }
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

            String sampleDir = mConfig.CobaltDir.replaceAll("\\*", sampleId);

            final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(sampleDir, sampleId);

            addCobaltSampleData(mConfig.RefGenVersion, amberGender, cobaltFilename, mChrRegionData);
        }
    }

    public static void main(final String[] args) throws ParseException
    {
        final Options options = new Options();
        NormalisationConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NormalisationFileBuilder normFileBuilder = new NormalisationFileBuilder(cmd);
        normFileBuilder.run();
    }

    private static CommandLine createCommandLine(final String[] args, final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
