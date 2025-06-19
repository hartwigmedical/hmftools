package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.score.NeoRnaData.NO_TPM_VALUE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Maps;

public class TpmMediansCache
{
    private final Map<String,Map<String,Double>> mTransCancerTpmMap;
    private final Set<String> mCancerTypes;

    private static final String ALL_TYPES = "All";
    public static final int CANCER_VALUE = 0;
    public static final int PAN_CANCER_VALUE = 1;

    public TpmMediansCache(final String filename)
    {
        mTransCancerTpmMap = Maps.newHashMap();
        mCancerTypes = Sets.newHashSet();

        if(filename != null)
            loadCohortFile(filename);
    }

    public boolean hasCancerType(final String cancerType) { return mCancerTypes.contains(cancerType); }

    public double[] getTranscriptTpm(final String transName, final String cancerType)
    {
        final double[] results = { NO_TPM_VALUE, NO_TPM_VALUE };
        final Map<String,Double> cancerTpmMap = mTransCancerTpmMap.get(transName);

        if(cancerTpmMap != null)
        {
            Double tpm = cancerTpmMap.get(ALL_TYPES);

            if(tpm != null)
                results[PAN_CANCER_VALUE] = tpm;

            if(cancerType != null && !cancerType.isEmpty() && cancerTpmMap.containsKey(cancerType))
                results[CANCER_VALUE] = cancerTpmMap.get(cancerType);
        }

        return results;
    }

    private static final String GENE_EXP_DELIM = ",";

    private void loadCohortFile(final String filename)
    {
        // GeneId,GeneName,TransName,Mesothelium,etc
        if(!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("invalid cohort TPM medians file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(header, GENE_EXP_DELIM);

            fieldsMap.keySet().stream()
                    .filter(x -> !x.equals(FLD_GENE_ID) && !x.equals(FLD_GENE_NAME) && !x.equals(FLD_TRANS_NAME))
                    .forEach(x -> mCancerTypes.add(x));

            int transNameIndex = fieldsMap.get(FLD_TRANS_NAME);

            String line = null;
            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(GENE_EXP_DELIM, -1);

                final String transName = items[transNameIndex];
                final Map<String,Double> cancerTpmMap = Maps.newHashMap();
                mTransCancerTpmMap.put(transName, cancerTpmMap);

                for(String cancerType : mCancerTypes)
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
        }
    }
}
