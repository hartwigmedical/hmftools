package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionData.RATE_COUNT;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionData.RATE_VALUE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_CHROMOSOME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;

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
    public final List<double[]> GeneRates;

    public GeneCohortData(final String geneId, final String geneName, final String chromosome)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;

        SampleIds = Lists.newArrayList();
        GeneRates = Lists.newArrayList();
    }

    public void addSampleData(final String sampleId, double geneRate)
    {
        int index = 0;
        while(index < GeneRates.size())
        {
            final double[] fpmData = GeneRates.get(index);

            if(Doubles.equal(geneRate, fpmData[RATE_VALUE]))
            {
                ++fpmData[RATE_COUNT];
                return;
            }

            if(geneRate < fpmData[RATE_VALUE])
                break;

            ++index;
        }

        SampleIds.add(index, sampleId);
        GeneRates.add(index, new double[] {geneRate, 1});
    }

    public static GeneCohortData fromCsv(final String[] items, final Map<String,Integer> fieldIndexMap)
    {
        return new GeneCohortData(
                items[fieldIndexMap.get(FLD_GENE_ID)],
                items[fieldIndexMap.get(FLD_GENE_NAME)],
                items[fieldIndexMap.get(FLD_CHROMOSOME)]);
    }

}
