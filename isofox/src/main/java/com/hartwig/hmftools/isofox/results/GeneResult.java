package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_SET_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.StringJoiner;

import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;

public class GeneResult
{
    public final com.hartwig.hmftools.common.gene.GeneData Gene;
    public final String CollectionId;
    public final int IntronicLength;
    public final int TransCount;

    private double mSplicedAlloc;
    private double mUnsplicedAlloc;
    private double mRawTpm;
    private double mAdjustedTpm;
    private double mFitResiduals;
    private double mLowMapQualsAllocation;

    public GeneResult(final GeneCollection geneCollection, final GeneReadData geneReadData)
    {
        Gene = geneReadData.GeneData;
        CollectionId = geneCollection.chrId();

        long exonicLength = geneReadData.calcExonicRegionLength();
        IntronicLength = (int)(Gene.length() - exonicLength);
        TransCount = geneReadData.getTranscripts().size();

        mFitResiduals = 0;
        mSplicedAlloc = 0;
        mRawTpm = 0;
        mAdjustedTpm = 0;
        mUnsplicedAlloc = 0;
        mLowMapQualsAllocation = 0;
    }

    public void setFitAllocation(double splicedAlloc, double unsplicedAlloc)
    {
        mSplicedAlloc = splicedAlloc;
        mUnsplicedAlloc = unsplicedAlloc;
    }

    public void setTPM(double raw, double adjusted)
    {
        mRawTpm = raw;
        mAdjustedTpm = adjusted;
    }

    public void applyTpmAdjustFactor(double factor) { mAdjustedTpm /= factor; }

    public void setFitResiduals(double residuals) { mFitResiduals = residuals; }
    public double getFitResiduals() { return mFitResiduals; }
    public double getSplicedAlloc() { return mSplicedAlloc; }
    public double getUnsplicedAlloc() { return mUnsplicedAlloc; }

    public void setLowMapQualsAllocation(double alloc) { mLowMapQualsAllocation = alloc; }

    public static String csvHeader()
    {
        return new StringJoiner(DELIMITER)
                .add(FLD_GENE_ID)
                .add(FLD_GENE_NAME)
                .add(FLD_CHROMOSOME)
                .add("GeneLength")
                .add("IntronicLength")
                .add("TranscriptCount")
                .add(FLD_GENE_SET_ID)
                .add(FLD_SPLICED_FRAGS)
                .add(FLD_UNSPLICED_FRAGS)
                .add(FLD_ADJ_TPM)
                .add("RawTPM")
                .add("FitResiduals")
                .add("LowMapQualFrags")
                .toString();
    }

    public String toCsv()
    {
        return new StringJoiner(DELIMITER)
                .add(Gene.GeneId)
                .add(Gene.GeneName)
                .add(Gene.Chromosome)
                .add(String.valueOf(Gene.length()))
                .add(String.valueOf(IntronicLength))
                .add(String.valueOf(TransCount))
                .add(CollectionId)
                .add(String.format("%.0f", mSplicedAlloc))
                .add(String.format("%.0f", mUnsplicedAlloc))
                .add(String.format("%6.3e", mAdjustedTpm))
                .add(String.format("%6.3e", mRawTpm))
                .add(String.format("%.1f", getFitResiduals()))
                .add(String.format("%.1f", mLowMapQualsAllocation))
                .toString();
    }
}
