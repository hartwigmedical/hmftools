package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;

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
        BufferedReader fileReader = createBufferedReader(cohortFile);

        String header = fileReader.readLine();
        String fileDelim = inferFileDelimiter(cohortFile);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
        int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
        int sampleCountIndex = fieldsIndexMap.get("SampleCount");

        String line = null;
        while((line = fileReader.readLine()) != null)
        {
            String[] items = line.split(fileDelim, -1);

            int sampleCount = Integer.parseInt(items[sampleCountIndex]);
            final String asjKey = formKey(items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

            mCohortFrequency.put(asjKey, sampleCount);
        }

        LOGGER.info(" Loaded {} alt-SJ sites from cohort file {}", mCohortFrequency.size(), cohortFile);
    }
}
