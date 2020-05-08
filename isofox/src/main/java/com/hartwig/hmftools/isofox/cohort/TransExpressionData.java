package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class TransExpressionData
{
    public final String GeneId;
    public final String GeneName;
    public final int TransId;
    public final String TransName;

    public final List<String> SampleIds;
    public final List<double[]> TpmValues;

    public static final int RATE_VALUE = 0;
    public static final int RATE_COUNT = 1;

    public TransExpressionData(final String geneId, final String geneName, final int transId, final String transName)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransId = transId;
        TransName = transName;

        SampleIds = Lists.newArrayList();
        TpmValues = Lists.newArrayList();
    }

    public void addSampleData(final String sampleId, double tpm)
    {
        int index = 0;
        while(index < TpmValues.size())
        {
            final double[] tpmData = TpmValues.get(index);

            if(Doubles.equal(tpm, tpmData[RATE_VALUE]))
            {
                ++tpmData[RATE_COUNT];
                return;
            }

            if(tpm < tpmData[RATE_VALUE])
                break;

            ++index;
        }

        SampleIds.add(index, sampleId);
        TpmValues.add(index, new double[] {tpm, 1});
    }


    public static TransExpressionData fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new TransExpressionData(
                items[fieldIndexMap.get(FLD_GENE_ID)],
                items[fieldIndexMap.get(FLD_GENE_NAME)],
                Integer.parseInt(items[fieldIndexMap.get(FLD_TRANS_ID)]),
                items[fieldIndexMap.get(FLD_TRANS_NAME)]);
    }

}
