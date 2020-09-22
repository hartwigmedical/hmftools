package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_SALMON;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.SOURCE_ISOFOX;

import java.util.Map;

public class ExpressionData
{
    public final String Source;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;

    private int mReadCount;
    private double mTPM;

    public ExpressionData(final String source, final String geneId, final String geneName, final String transName, int readCount, double tpm)
    {
        Source = source;
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        mReadCount = readCount;
        mTPM = tpm;
    }

    public void addCounts(double tpm, int readCount)
    {
        mTPM += tpm;
        mReadCount += readCount;
    }

    public void setTpm(double tpm) { mTPM = tpm; }
    public double tpm() { return mTPM; }
    public int readCount() { return mReadCount; }

    public static ExpressionData fromIsofoxTranscript(final String data, int geneIdIndex, int geneNameIndex, int transIndex, int readCountIndex, int tpmIndex)
    {
        final String[] items = data.split(",");

        return new ExpressionData(
                SOURCE_ISOFOX, items[geneIdIndex], items[geneNameIndex], items[transIndex],
                Integer.parseInt(items[readCountIndex]), Double.parseDouble(items[tpmIndex]));
    }

    public static ExpressionData fromIsofoxGene(final String data, int geneIdIndex, int geneNameIndex, int readCountIndex, int tpmIndex)
    {
        final String[] items = data.split(",");

        return new ExpressionData(
                SOURCE_ISOFOX, items[geneIdIndex], items[geneNameIndex], "",
                Integer.parseInt(items[readCountIndex]), Double.parseDouble(items[tpmIndex]));
    }

    public static ExpressionData fromSalmon(final String data, final Map<String,String[]> geneTransMap)
    {
        // Name    Length  EffectiveLength TPM     NumReads
        //ENST00000415118.1       8       9.000   0.000000        0.000

        final String[] items = data.split("\t");

        if(items.length != 5)
            return null;

        String transName = items[0].replaceAll("\\.[0-9]*",""); // strip off trans index

        final String[] geneData = geneTransMap.get(transName);

        if(geneData == null)
            return null;

        final String geneId = geneData != null ? geneData[0] : "";
        final String geneName = geneData != null ? geneData[1] : "";

        double tpm = Double.parseDouble(items[3]);
        int readCount = (int)Double.parseDouble(items[4]);

        return new ExpressionData(EXT_SOURCE_SALMON, geneId, geneName, transName, readCount, tpm);
    }
}
