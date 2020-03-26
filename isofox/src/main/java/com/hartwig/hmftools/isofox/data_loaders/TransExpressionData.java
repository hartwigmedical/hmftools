package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.Map;

public class TransExpressionData
{
    public final String GeneId;
    public final String GeneName;
    public final int TransId;
    public final String TransName;
    public final int EffectiveLength;
    public final double FitAllocation;

    private double mTransPerM;
    private double mFragsPerKb;

    public TransExpressionData(final String geneId, final String geneName, final int transId, final String transName,
            final int effectiveLength, final double fitAllocation)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransId = transId;
        TransName = transName;
        EffectiveLength = effectiveLength;
        FitAllocation = fitAllocation;

        mTransPerM = 0;
        mFragsPerKb = effectiveLength > 0 ? fitAllocation / (effectiveLength / 1000.0) : 0;
    }

    public double fragsPerKb() { return mFragsPerKb; }

    public double transPerM() { return mTransPerM; }
    public void setTransPerM(double tpm) { mTransPerM = tpm; }

    public static TransExpressionData fromCsv(final String input, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = input.split(DELIMITER);

        return new TransExpressionData(
                items[fieldIndexMap.get("GeneId")],
                items[fieldIndexMap.get("GeneName")],
                Integer.parseInt(items[fieldIndexMap.get("TransId")]),
                items[fieldIndexMap.get("TransName")],
                Integer.parseInt(items[fieldIndexMap.get("EffectiveLength")]),
                Double.parseDouble(items[fieldIndexMap.get("FitAllocation")]));
    }

}
