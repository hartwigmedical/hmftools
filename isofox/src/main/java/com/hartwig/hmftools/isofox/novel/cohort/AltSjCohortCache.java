package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.FileReader;
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
            BufferedReader fileReader = new BufferedReader(new FileReader(cohortFile));

            String line = fileReader.readLine();
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            int geneIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int sampleCountIndex = fieldsIndexMap.get("SampleCount");

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneId = items[geneIndex];

                // if(!mConfig.processGeneId(geneId))
                //    continue;

                int sampleCount = Integer.parseInt(items[sampleCountIndex]);
                final String asjKey = formKey(items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

                mCohortFrequency.put(asjKey, sampleCount);
            }

            ISF_LOGGER.info("loaded alt-SJ cohort file({}) with {} sites", cohortFile, mCohortFrequency.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load alt-SJ cohort file({}): {}", cohortFile, e.toString());
            return;
        }
    }

}
