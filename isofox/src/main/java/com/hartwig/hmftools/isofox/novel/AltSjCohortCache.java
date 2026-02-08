package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;

import com.google.common.collect.Maps;

public class AltSjCohortCache
{
    private final Map<String,Integer> mCohortFrequency;

    public AltSjCohortCache(final String cohortFile)
    {
        mCohortFrequency = Maps.newHashMap();
        loadAltSjCohortFile(cohortFile);
    }

    public int getCohortFrequency(final String asjKey)
    {
        Integer cohortFrequency = mCohortFrequency.get(asjKey);
        return cohortFrequency != null ? cohortFrequency : 1;
    }

    public static String formKey(final String chromosome, int posStart, int posEnd)
    {
        return String.format("%s;%s;%s", chromosome, posStart, posEnd);
    }

    private void loadAltSjCohortFile(final String cohortFile)
    {
        if(cohortFile == null)
            return;

        if(!Files.exists(Paths.get(cohortFile)))
        {
            ISF_LOGGER.error("missing alt-SJ cohort file");
            return;
        }

        try
        {
            BufferedReader fileReader = createBufferedReader(cohortFile);

            String line = fileReader.readLine();
            String fileDelim = inferFileDelimiter(cohortFile);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, fileDelim);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int sampleCountIndex = fieldsIndexMap.get("SampleCount");

            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim, -1);

                int sampleCount = Integer.parseInt(values[sampleCountIndex]);
                final String asjKey = formKey(values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex]));

                mCohortFrequency.put(asjKey, sampleCount);
            }

            ISF_LOGGER.debug("loaded {} alt-SJ sites from cohort file({})", mCohortFrequency.size(), cohortFile);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load alt-SJ cohort file({}): {}", cohortFile, e.toString());
        }
    }
}
