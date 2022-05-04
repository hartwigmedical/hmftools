package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
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

            if(cancerType != null && !cancerType.isEmpty() && cancerTpmMap.containsKey(cancerType))
            {
                results[CANCER_VALUE] = cancerTpmMap.get(cancerType);
            }
        }

        return results;
    }

    private static final String GENE_EXP_DELIM = ",";
    private static final String FLD_GENE_ID = "GeneId";
    private static final String FLD_GENE_NAME = "GeneName";
    private static final String FLD_TRANS_NAME = "TransName";


    private void loadCohortFile(final String filename)
    {
        // GeneId,GeneName,TransName,Mesothelium,etc
        if(!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("invalid cohort TPM file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
                return;

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(line, GENE_EXP_DELIM);

            final List<String> cancerTypes = fieldsMap.keySet().stream()
                    .filter(x -> !x.equals(FLD_GENE_ID) && !x.equals(FLD_GENE_NAME) && !x.equals(FLD_TRANS_NAME))
                    .collect(Collectors.toList());

            int transNameIndex = fieldsMap.get(FLD_TRANS_NAME);

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(GENE_EXP_DELIM, -1);

                final String transName = items[transNameIndex];
                final Map<String,Double> cancerTpmMap = Maps.newHashMap();
                mTransCancerTpmMap.put(transName, cancerTpmMap);

                for(String cancerType : cancerTypes)
                {
                    double tpm = Double.parseDouble(items[fieldsMap.get(cancerType)]);
                    cancerTpmMap.put(cancerType, tpm);
                }
            }

            NE_LOGGER.info("loaded {} cohort TPM median data", mTransCancerTpmMap.size());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load cohort TPM median data file({}): {}", filename, e.toString());
            return;
        }

    }
}
