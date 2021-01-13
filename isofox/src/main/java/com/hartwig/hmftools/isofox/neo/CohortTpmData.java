package com.hartwig.hmftools.isofox.neo;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class CohortTpmData
{
    private Map<String,Map<String,Double>> mTransCancerTpmMap;

    private static final String ALL_TYPES = "All";
    public static final int CANCER_VALUE = 0;
    public static final int COHORT_VALUE = 1;

    public CohortTpmData(final String filename)
    {
        mTransCancerTpmMap = Maps.newHashMap();

        if(filename != null)
            loadCohortFile(filename);
    }

    public double[] getTranscriptTpm(final String transName, final String cancerType)
    {
        final double[] results = { 0, 0 };
        final Map<String,Double> cancerTpmMap = mTransCancerTpmMap.get(transName);

        if(cancerTpmMap != null)
        {
            Double tpm = cancerTpmMap.get(ALL_TYPES);
            results[COHORT_VALUE] = tpm != null ? tpm : 0;

            if(cancerType != null && cancerTpmMap.containsKey(cancerType))
            {
                results[CANCER_VALUE] = cancerTpmMap.get(cancerType);
            }
        }

        return results;
    }

    private void loadCohortFile(final String filename)
    {
        // GeneId,GeneName,TransName,Mesothelium,etc
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid cohort TPM file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
                return;

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(line, DELIMITER);

            final List<String> cancerTypes = fieldsMap.keySet().stream()
                    .filter(x -> !x.equals(FLD_GENE_ID) && !x.equals(FLD_GENE_NAME) && !x.equals(FLD_TRANS_NAME))
                    .collect(Collectors.toList());

            int transNameIndex = fieldsMap.get(FLD_TRANS_NAME);

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String transName = items[transNameIndex];
                final Map<String,Double> cancerTpmMap = Maps.newHashMap();
                mTransCancerTpmMap.put(transName, cancerTpmMap);

                for(String cancerType : cancerTypes)
                {
                    double tpm = Double.parseDouble(items[fieldsMap.get(cancerType)]);
                    cancerTpmMap.put(cancerType, tpm);
                }
            }

            ISF_LOGGER.info("loaded {} cohort TPM median data", mTransCancerTpmMap.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load cohort TPM median data file({}): {}", filename, e.toString());
            return;
        }

    }
}
