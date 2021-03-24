package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.utils.Matrix;

public class SomaticSigs
{
    private final List<String> mSignatureNames;
    private final Matrix mSignatures;
    private final LeastSquaresFit mLeastSquaresFitter;

    public SomaticSigs(final String signaturesFile)
    {
        mSignatureNames = Lists.newArrayList();

        if(signaturesFile != null && !signaturesFile.isEmpty() && Files.exists(Paths.get(signaturesFile)))
        {
            mSignatures = loadMatrixDataFile(signaturesFile, mSignatureNames);
            mLeastSquaresFitter = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);
        }
        else
        {
            mSignatures = null;
            mLeastSquaresFitter = null;
        }
    }

    public boolean hasValidData() { return mSignatures != null && !mSignatureNames.isEmpty(); }

    public String getSigName(int index) { return index < mSignatureNames.size() ? mSignatureNames.get(index) : ""; }

    public final double[] fitSampleCounts(final double[] sampleCounts)
    {
        double sampleTotal = sumVector(sampleCounts);

        if(sampleTotal == 0)
            return null;

        mLeastSquaresFitter.initialise(mSignatures.getData(), sampleCounts);
        mLeastSquaresFitter.solve();
        return mLeastSquaresFitter.getContribs();
    }
}
