package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.positionKey;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;


public class UmrCohortFrequency
{
    private final Map<String,Map<String,Integer>> mCohortFrequency; // keyed by chromosome and then pos-key

    public static final String UMR_COHORT_FREQUENCY_FILE = "umr_cohort_freq_file";

    public UmrCohortFrequency(final String cohortFile)
    {
        mCohortFrequency = Maps.newHashMap();
        loadFile(cohortFile);
    }

    public int getCohortFrequency(final String chromosome, final int scSide, final int exonBoundary)
    {
        return getCohortFrequency(chromosome, positionKey(scSide, exonBoundary));
    }

    public int getCohortFrequency(final String chromosome, final String positionKey)
    {
        Map<String,Integer> locationMap = mCohortFrequency.get(chromosome);

        if(locationMap == null)
            return 0;

        Integer frequency = locationMap.get(positionKey);
        return frequency != null ? frequency : 0;
    }

    private void loadFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int scSideIndex = fieldsIndexMap.get("SoftClipSide");
            int exonBoundaryIndex = fieldsIndexMap.get("ExonBoundary");
            int cohortFreqIndex = fieldsIndexMap.get("CohortFreq");

            String line = "";
            int count = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIMITER, -1);

                String chromosome = values[chrIndex];

                int scSide = Integer.parseInt(values[scSideIndex]);
                int exonBoundary = Integer.parseInt(values[exonBoundaryIndex]);
                int cohortFreq = Integer.parseInt(values[cohortFreqIndex]);

                Map<String,Integer> locationMap = mCohortFrequency.get(chromosome);

                if(locationMap == null)
                {
                    locationMap = Maps.newHashMap();
                    mCohortFrequency.put(chromosome, locationMap);
                }

                locationMap.put(positionKey(scSide, exonBoundary), cohortFreq);
                ++count;
            }

            ISF_LOGGER.info("loaded {} unmapped-read cohort frequency records from {}", count, filename);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load unmapped-read cohort frequency file({}): {}", filename.toString(), e.toString());
        }
    }
}
