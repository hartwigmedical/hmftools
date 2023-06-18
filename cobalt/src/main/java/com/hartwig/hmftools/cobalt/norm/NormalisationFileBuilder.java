package com.hartwig.hmftools.cobalt.utils;

import static java.lang.Math.floor;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.utils.CobaltDataLoader.addCobaltSampleData;
import static com.hartwig.hmftools.cobalt.utils.NormConstants.MIN_ENRICHMENT_RATIO;
import static com.hartwig.hmftools.cobalt.utils.NormConstants.REGION_SIZE;
import static com.hartwig.hmftools.cobalt.utils.Normaliser.calcRelativeEnrichment;
import static com.hartwig.hmftools.cobalt.utils.Normaliser.calcSampleAdjustedRatios;
import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.genome.gc.GCBucket.calcGcBucket;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberGender;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.purple.Gender;

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

        // load per-sample Cobalt ratios
        loadSampleCobaltData(sampleGenders);

        // calculate per-sample normalised tumor GC ratios
        calcSampleAdjustedRatios(mConfig.SampleIds, mChrRegionData);

        // calculate a final relative panel enichment ratio for each region
        calcRelativeEnrichment(mChrRegionData, MIN_ENRICHMENT_RATIO);

        writeNormalisationFile();
        writeDetailedFile();

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

                addRegion(regions, startPosition);

                while(endPosition < namedBed.end())
                {
                    startPosition += REGION_SIZE;
                    addRegion(regions, startPosition);
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

    private static void addRegion(final List<RegionData> regions, int position)
    {
        RegionData prevRegion = !regions.isEmpty() ? regions.get(regions.size() - 1) : null;
        if(prevRegion != null && prevRegion.Position == position)
            return;

        regions.add(new RegionData(position));
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

    private void writeNormalisationFile()
    {
        CB_LOGGER.info("writing normalisation file: {}", mConfig.OutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputFile, false);

            writer.write("chromosome\tposition\trelativeEnrichment");
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                if(!mChrRegionData.containsKey(chrStr))
                    continue;

                for(RegionData regionData : mChrRegionData.get(chrStr))
                {
                    writer.write(format("%s\t%d\t%.4f",
                            chrStr, regionData.Position,
                            regionData.relativeEnrichment() > 0 ? regionData.relativeEnrichment() : Double.NaN));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write normalisation file: {}", e.toString());
            System.exit(1);
        }
    }

    private void writeDetailedFile()
    {
        if(mConfig.DetailedFile == null)
            return;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.DetailedFile, false);

            writer.write("SampleId\tChromosome\tPosition\tGcBucket\tMappability\tGcRatioPanel\tReadCount\tAdjGcRatio");
            writer.newLine();

            for(Map.Entry<String, List<RegionData>> entry : mChrRegionData.entrySet())
            {
                String chromosome = entry.getKey();

                List<RegionData> regions = entry.getValue();

                for(RegionData regionData : regions)
                {
                    for(int i = 0; i < mConfig.SampleIds.size(); ++i)
                    {
                        String sampleId = mConfig.SampleIds.get(i);
                        SampleRegionData sampleRegionData = regionData.getSampleData(i);

                        writer.write(format("%s\t%s\t%d\t%d\t%.3f\t%.3f\t%d\t%.3f",
                                sampleId, chromosome, regionData.Position, regionData.gcBucket(), regionData.mappability(),
                                sampleRegionData.GcRatioPanel, sampleRegionData.ReadCount, sampleRegionData.adjustedGcRatio()));
                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write normalisation file: {}", e.toString());
            System.exit(1);
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
