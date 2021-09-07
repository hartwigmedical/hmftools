package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_RSEM;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_SALMON;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.SOURCE_ISOFOX;

import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;

public class ExpressionData
{
    public final String Source;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final int EffectiveLength;

    private double mFittedFragmentCount;
    private double mRawFragmentCount;
    private int mReadCount;
    private double mTPM;

    public final int SplicedFragments;
    public final int UnsplicedFragments;
    public final double LowMapQualFrags;

    public ExpressionData(
            final String source, final String geneId, final String geneName, final String transName,
            double fittedFrags, double rawFrags, int readCount, double tpm, int effectiveLength,
            int splicedFragments, int unsplicedFragments, double lowMapQualFrags)
    {
        Source = source;
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        EffectiveLength = effectiveLength;
        SplicedFragments = splicedFragments;
        UnsplicedFragments = unsplicedFragments;
        LowMapQualFrags = lowMapQualFrags;

        mFittedFragmentCount = fittedFrags;
        mRawFragmentCount = rawFrags;
        mReadCount = readCount;
        mTPM = tpm;
    }

    public void addCounts(double tpm, double fittedFrags, double rawFrags, int reads)
    {
        mTPM += tpm;
        mReadCount += reads;
        mRawFragmentCount += rawFrags;
        mFittedFragmentCount += fittedFrags;
    }

    public double tpm() { return mTPM; }
    public int readCount() { return mReadCount; }
    public double fittedFragments() { return mFittedFragmentCount; }
    public double rawFragment() { return mRawFragmentCount; }

    public static ExpressionData fromIsofoxTranscript(
            final String data, int geneIdIndex, int geneNameIndex, int transIndex,
            int fittedFragIndex, int rawFragsIndex, int tpmIndex, int effectiveLengthIndex, int lowQualIndex)
    {
        final String[] items = data.split(",");

        return new ExpressionData(
                SOURCE_ISOFOX, items[geneIdIndex], items[geneNameIndex], items[transIndex],
                Double.parseDouble(items[fittedFragIndex]), Double.parseDouble(items[rawFragsIndex]),
                0, Double.parseDouble(items[tpmIndex]), Integer.parseInt(items[effectiveLengthIndex]),
                0, 0, lowQualIndex >= 0 ? Double.parseDouble(items[lowQualIndex]) : 0);
    }

    public static ExpressionData fromIsofoxGene(
            final String data, int geneIdIndex, int geneNameIndex, int tpmIndex, int splicedIndex, int unsplicedIndex, int lowQualIndex)
    {
        final String[] items = data.split(",");

        return new ExpressionData(
                SOURCE_ISOFOX, items[geneIdIndex], items[geneNameIndex], "",
                0, 0, 0, Double.parseDouble(items[tpmIndex]), 0,
                Integer.parseInt(items[splicedIndex]), Integer.parseInt(items[unsplicedIndex]),
                lowQualIndex >= 0 ? Double.parseDouble(items[lowQualIndex]) : 0);
    }

    public static String getExternalSourceFilename(final String source, final String sampleId, boolean transScope)
    {
        if(source.equals(EXT_SOURCE_SALMON))
        {
            if(transScope)
                return String.format("%s.salmon.tsv", sampleId);
            else
                return null;
        }
        else
        {
            if(transScope)
                return String.format("%s.rsem.trans_data.tsv", sampleId);
            else
            return String.format("%s.rsem.gene_data.tsv", sampleId);
        }
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
        int effectiveLength = (int)Double.parseDouble(items[2]);

        return new ExpressionData(
                EXT_SOURCE_SALMON, geneId, geneName, transName, 0, 0, readCount, tpm,
                effectiveLength, 0, 0, 0);
    }

    public static ExpressionData fromRsemTranscript(final String data, final Map<String,String[]> geneTransMap)
    {
        // transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
        //  ENST00000373020.4       ENSG00000000003.10      2206    2040.02 167.87  1.50    2.78    81.71

        final String[] items = data.split("\t");

        if(items.length != 8)
            return null;

        String transName = items[0].replaceAll("\\.[0-9]*",""); // strip off index
        String geneId = items[1].replaceAll("\\.[0-9]*","");

        final String[] geneData = geneTransMap.get(transName);

        if(geneData == null)
            return null;

        final String geneName = geneData != null ? geneData[1] : "";

        double tpm = Double.parseDouble(items[5]);
        int readCount = 0;
        int effectiveLength = (int)Double.parseDouble(items[3]);

        return new ExpressionData(
                EXT_SOURCE_RSEM, geneId, geneName, transName, 0, 0, readCount, tpm,
                effectiveLength, 0, 0, 0);
    }

    public static ExpressionData fromRsemGene(final String data, final EnsemblDataCache geneTransCache)
    {
        // gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
        // ENSG00000000003.10      ENST00000373020.4,ENST00000494424.1,ENST00000496771.1   1963.14 1797.16 181.00  1.84    3.40
        // ENSG00000000005.5       ENST00000373031.4,ENST00000485971.1     940.50  774.56  0.00    0.00    0.00

        final String[] items = data.split("\t");

        if(items.length != 7)
            return null;

        String geneId = items[0].replaceAll("\\.[0-9]*","");

        final GeneData geneData = geneTransCache.getGeneDataById(geneId);

        if(geneData == null)
            return null;

        final String geneName = geneData.GeneName;

        double tpm = Double.parseDouble(items[5]);
        int readCount = 0;
        int effectiveLength = (int)Double.parseDouble(items[3]);

        return new ExpressionData(
                EXT_SOURCE_RSEM, geneId, geneName, "", 0, 0, readCount, tpm,
                effectiveLength, 0, 0, 0);
    }
}
