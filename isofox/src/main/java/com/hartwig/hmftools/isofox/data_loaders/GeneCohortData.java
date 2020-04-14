package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionData.TPM_COUNT;
import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionData.TPM_VALUE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_CHROMOSOME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class GeneCohortData
{
    public final String GeneId;
    public final String GeneName;
    public final String Chromosome;

    public final List<String> SampleIds;
    public final List<double[]> FragsPerMillionValues;

    public GeneCohortData(final String geneId, final String geneName, final String chromosome)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;

        SampleIds = Lists.newArrayList();
        FragsPerMillionValues = Lists.newArrayList();
    }

    public void addSampleData(final String sampleId, double fragsPerMillion)
    {
        int index = 0;
        while(index < FragsPerMillionValues.size())
        {
            final double[] fpmData = FragsPerMillionValues.get(index);

            if(Doubles.equal(fragsPerMillion, fpmData[TPM_VALUE]))
            {
                ++fpmData[TPM_COUNT];
                return;
            }

            if(fragsPerMillion < fpmData[TPM_VALUE])
                break;

            ++index;
        }

        SampleIds.add(index, sampleId);
        FragsPerMillionValues.add(index, new double[] {fragsPerMillion, 1});
    }

    public static GeneCohortData fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new GeneCohortData(
                items[fieldIndexMap.get(FLD_GENE_ID)],
                items[fieldIndexMap.get(FLD_GENE_NAME)],
                items[fieldIndexMap.get(FLD_CHROMOSOME)]);
    }

}
