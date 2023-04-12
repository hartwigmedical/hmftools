package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;

import java.util.List;

import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigResiduals;

public class CnPurityCalculator
{
    private double mFitCoefficient;
    private double mFitIntercept;
    private double mEstimatedPurity;
    private double mResiduals;
    private boolean mValid;

    public CnPurityCalculator()
    {
        mFitCoefficient = 0;
        mFitIntercept = 0;
        mEstimatedPurity = 0;
        mResiduals = 0;
        mValid = false;
    }

    public double fitCoefficient() { return mFitCoefficient; }
    public double fitIntercept() { return mFitIntercept; }
    public double estimatedPurity() { return mEstimatedPurity; }
    public double residuals() { return mResiduals; }
    public boolean valid() { return mValid; }

    public void calculatePurity(final List<CopyNumberGcData> copyNumberSegments, final double samplePloidy)
    {
        double gcRatioCountTotal = 0;
        double gcRatioMedianCountTotal = 0;
        double copyNumberCountTotal = 0;

        for(CopyNumberGcData cnSegment : copyNumberSegments)
        {
            gcRatioCountTotal += cnSegment.count();
            copyNumberCountTotal += cnSegment.count() * cnSegment.CopyNumber;
            gcRatioMedianCountTotal += cnSegment.count() * cnSegment.median();
        }

        double medianOffset = gcRatioMedianCountTotal / gcRatioCountTotal;
        double copyNumberOffset = copyNumberCountTotal / gcRatioCountTotal;
        int segmentCount = copyNumberSegments.size();

        double[][] adjustedCopyNumber = new double[segmentCount][1];
        double[] adjustedGcRatioMedians = new double[segmentCount];

        for(int i = 0; i < segmentCount; ++i)
        {
            CopyNumberGcData cnSegment = copyNumberSegments.get(i);

            double sqrtCount = sqrt(cnSegment.count());
            adjustedCopyNumber[i][0] = (cnSegment.CopyNumber - copyNumberOffset) * sqrtCount;

            double adjustMedian = (cnSegment.median() - medianOffset) * sqrtCount;
            adjustedGcRatioMedians[i] = adjustMedian;
        }

        LeastSquaresFit lsqFit = new LeastSquaresFit(segmentCount, 1);

        lsqFit.initialise(adjustedCopyNumber, adjustedGcRatioMedians);
        lsqFit.solve();

        mFitCoefficient = lsqFit.getContribs()[0];
        mFitIntercept = medianOffset - copyNumberOffset * mFitCoefficient;

        double fittedDiffTotal = 0;
        double fittedDiffTotalAbs = 0;
        double gcRatioTotal = 0;

        for(int i = 0; i < segmentCount; ++i)
        {
            CopyNumberGcData cnSegment = copyNumberSegments.get(i);
            double fittedGcRatio = mFitIntercept + cnSegment.CopyNumber * mFitCoefficient;
            gcRatioTotal += cnSegment.median();
            double gcRatioFitDiff = cnSegment.median() - fittedGcRatio;
            fittedDiffTotal += gcRatioFitDiff; // sum of LSQ fit
            fittedDiffTotalAbs += abs(gcRatioFitDiff);
        }

        mResiduals = gcRatioTotal > 0 ? fittedDiffTotalAbs / gcRatioTotal : 0;

        mEstimatedPurity = 2 * mFitCoefficient / (1 + (2 - samplePloidy) * mFitCoefficient);

        mValid = true;
    }

}
