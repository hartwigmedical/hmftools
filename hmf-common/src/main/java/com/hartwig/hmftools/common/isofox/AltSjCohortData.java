package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.RnaCommon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class AltSjCohortData
{
    private static final Logger LOGGER = LogManager.getLogger(AltSjCohortData.class);

    private final Map<String, Integer> mCohortFrequency = Maps.newHashMap();

    public AltSjCohortData(final String cohortFile) throws IOException
    {
        loadAltSjCohortFile(cohortFile);
    }

    public int getCohortFrequency(final String asjKey)
    {
        Integer cohortFrequency = mCohortFrequency.get(asjKey);
        return cohortFrequency != null ? cohortFrequency : 1;
    }

    private void loadAltSjCohortFile(final String cohortFile) throws IOException
    {
        BufferedReader fileReader = new BufferedReader(new FileReader(cohortFile));

        String line = fileReader.readLine();
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, RnaCommon.DELIMITER);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
        int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
        int sampleCountIndex = fieldsIndexMap.get("SampleCount");

        while((line = fileReader.readLine()) != null)
        {
            String[] items = line.split(RnaCommon.DELIMITER, -1);

            int sampleCount = Integer.parseInt(items[sampleCountIndex]);
            final String asjKey = formKey(items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

            mCohortFrequency.put(asjKey, sampleCount);
        }

        LOGGER.info(" Loaded {} alt-SJ sites from cohort file {}", mCohortFrequency.size(), cohortFile);
    }
}
