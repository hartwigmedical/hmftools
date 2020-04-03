package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_EFFECTIVE_LENGTH;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FIT_ALLOCATION;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class TransExpressionData
{
    public final String GeneId;
    public final String GeneName;
    public final int TransId;
    public final String TransName;

    public final List<String> SampleIds;
    public final List<Double> FitAllocations;
    public final List<Double> TpmValues;
    public final List<Integer> EffectiveLengths;

    public TransExpressionData(final String geneId, final String geneName, final int transId, final String transName)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransId = transId;
        TransName = transName;

        SampleIds = Lists.newArrayList();
        FitAllocations = Lists.newArrayList();
        TpmValues = Lists.newArrayList();
        EffectiveLengths = Lists.newArrayList();
    }

    public void addSampleData(final String sampleId, double fitAllocation, double tpm, int effectiveLength)
    {
        int index = 0;
        while(index < TpmValues.size())
        {
            if(tpm < TpmValues.get(index))
                break;

            ++index;
        }

        SampleIds.add(index, sampleId);
        FitAllocations.add(index, fitAllocation);
        TpmValues.add(index, tpm);
        EffectiveLengths.add(index, effectiveLength);
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
