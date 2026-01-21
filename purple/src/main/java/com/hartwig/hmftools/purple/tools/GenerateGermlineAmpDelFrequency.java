package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.FLD_FREQUENCY;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.FLD_TYPE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;

import org.jetbrains.annotations.NotNull;

public class GenerateGermlineAmpDelFrequency
{
    private final String mCohortFrequencyFile;
    private final String mCohortDeletionsFile;
    private final int mMinSampleCount;
    private final List<String> mSampleIds;

    private Map<String,List<AmpDelRegionFrequency>> mChrRegionMap;

    public static final String COHORT_AMP_DEL_FREQ_FILE = "germline_amp_del_freq_file";

    private static final String MIN_SAMPLE_COUNT = "min_samples";
    private static final String COHORT_AMP_DEL_FILE = "cohort_germline_amp_del_file";
    private static final int DEFAULT_MIN_SAMPLES = 3;

    public GenerateGermlineAmpDelFrequency(final ConfigBuilder configBuilder)
    {
        mCohortDeletionsFile = configBuilder.getValue(COHORT_AMP_DEL_FILE);
        mCohortFrequencyFile = configBuilder.getValue(COHORT_AMP_DEL_FREQ_FILE);
        mMinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);

        mChrRegionMap = Maps.newHashMap();

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            mSampleIds = loadSampleIdsFile(configBuilder);
            PPL_LOGGER.info("loaded {} reference samples", mSampleIds.size());
        }
        else
        {
            mSampleIds = Lists.newArrayList();
        }
    }

    public void buildCohortFrequencies()
    {
        if(mCohortDeletionsFile == null || mCohortFrequencyFile == null)
        {
            PPL_LOGGER.error("missing germline deletion file or frequency output file not configured");
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mCohortDeletionsFile));
            String header = fileReader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            int sampleIdIndex = fieldsIndexMap.get("SampleId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int regionStartIndex = fieldsIndexMap.get("RegionStart");
            int regionEndIndex = fieldsIndexMap.get("RegionEnd");
            int germlineStatusIndex = fieldsIndexMap.get("GermlineStatus");

            String line = "";
            int count = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                String sampleId = values[sampleIdIndex];

                if(!mSampleIds.isEmpty() && !mSampleIds.contains(sampleId))
                    continue;

                String chromosome = values[chromosomeIndex];
                int regionStart = Integer.parseInt(values[regionStartIndex]);
                int regionEnd = Integer.parseInt(values[regionEndIndex]);
                GermlineStatus germlineStatus = GermlineStatus.valueOf(values[germlineStatusIndex]);

                addEvent(chromosome, regionStart, regionEnd, germlineStatus);

                ++count;

                if(count > 0 && (count % 100000) == 0)
                {
                    PPL_LOGGER.debug("processed {} cohort germline deletions", count);
                }
            }

            PPL_LOGGER.info("processed {} cohort germline deletions", count);
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to read cohort germline deletions file: {}", e.toString());
        }

        writeCohortFrequencies();
    }

    public void addEvent(final String chromsome, int regionStart, int regionEnd, GermlineStatus germlineStatus)
    {
        List<AmpDelRegionFrequency> regions = mChrRegionMap.get(chromsome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrRegionMap.put(chromsome, regions);
        }
        AmpDelRegionFrequency.EventType type = germlineStatus == GermlineStatus.AMPLIFICATION ? AmpDelRegionFrequency.EventType.AMP : AmpDelRegionFrequency.EventType.DEL;

        int index = 0;
        while(index < regions.size())
        {
            AmpDelRegionFrequency region = regions.get(index);

            if(region.Region.start() == regionStart && region.Region.end() == regionEnd && region.Type.equals(type))
            {
                ++region.Frequency;
                return;
            }

            if(regionStart < region.Region.start() || (regionStart == region.Region.start() && regionEnd < region.Region.end()))
            {
                break;
            }

            ++index;
        }
        regions.add(index, new AmpDelRegionFrequency(new BaseRegion(regionStart, regionEnd), type, 1));
    }

    private void writeCohortFrequencies()
    {
        PPL_LOGGER.info("writing cohort germline frequencies file: {}", mCohortFrequencyFile);

        try
        {

            BufferedWriter writer = createBufferedWriter(mCohortFrequencyFile, false);

            StringJoiner sj = new StringJoiner(CSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_REGION_START).add(FLD_REGION_END).add(FLD_TYPE).add(FLD_FREQUENCY);
            writer.write(sj.toString());
            writer.newLine();

            for(Map.Entry<String,List<AmpDelRegionFrequency>> entry : mChrRegionMap.entrySet())
            {
                String chromosome = entry.getKey();

                for(AmpDelRegionFrequency region : entry.getValue())
                {
                    if(region.Frequency < mMinSampleCount)
                        continue;

                    writer.write(String.format("%s,%d,%d,%s,%d", chromosome, region.Region.start(), region.Region.end(), region.Type, region.Frequency));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write cohort germline deletion frequency file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(COHORT_AMP_DEL_FILE, true, "Input germline cohort amps/dels file");
        configBuilder.addPath(COHORT_AMP_DEL_FREQ_FILE, false, "Output cohort germline amp/del frequency file");
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min sample frequency to write a region", DEFAULT_MIN_SAMPLES);
        configBuilder.addConfigItem(SAMPLE_ID_FILE, true, "Reference de-duped sample IDs");
        addLoggingOptions(configBuilder);

        setLogLevel(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateGermlineAmpDelFrequency germlineDeletionFrequency = new GenerateGermlineAmpDelFrequency(configBuilder);
        germlineDeletionFrequency.buildCohortFrequencies();
    }
}
