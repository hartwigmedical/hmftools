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

        // compute average GC ratio median and average copy number
        int segmentCount = 0;

        for(CopyNumberGcData cnSegment : copyNumberSegments)
        {
            if(!cnSegment.IsValid)
                continue;

            ++segmentCount;
            gcRatioCountTotal += cnSegment.count();
            copyNumberCountTotal += cnSegment.count() * cnSegment.CopyNumberLevel;
            gcRatioMedianCountTotal += cnSegment.count() * cnSegment.median();
        }

        double avgGcRatio = gcRatioMedianCountTotal / gcRatioCountTotal;
        double avgCopyNumber = copyNumberCountTotal / gcRatioCountTotal;

        double[][] adjustedCopyNumber = new double[segmentCount][1];
        double[] adjustedGcRatioMedians = new double[segmentCount];

        // weight the data by the number of GC ratio points in each Cn segment
        int segmentIndex = 0;
        for(CopyNumberGcData cnSegment : copyNumberSegments)
        {
            if(!cnSegment.IsValid)
                continue;

            double sqrtCount = sqrt(cnSegment.count());
            adjustedCopyNumber[segmentIndex][0] = (cnSegment.CopyNumberLevel - avgCopyNumber) * sqrtCount;

            double adjustMedian = (cnSegment.median() - avgGcRatio) * sqrtCount;
            adjustedGcRatioMedians[segmentIndex] = adjustMedian;

            ++segmentIndex;
        }

        LeastSquaresFit lsqFit = new LeastSquaresFit(segmentCount, 1);

        lsqFit.initialise(adjustedCopyNumber, adjustedGcRatioMedians);
        lsqFit.solve();

        mFitCoefficient = lsqFit.getContribs()[0];
        mFitIntercept = avgGcRatio - avgCopyNumber * mFitCoefficient;

        double fittedDiffWeightedTotalAbs = 0;
        double gcRatioWeightedTotal = 0;

        for(CopyNumberGcData cnSegment : copyNumberSegments)
        {
            if(!cnSegment.IsValid)
                continue;

            double fittedGcRatio = mFitIntercept + cnSegment.CopyNumberLevel * mFitCoefficient;
            double gcRatioFitDiff = cnSegment.median() - fittedGcRatio;

            double sqrtCount = sqrt(cnSegment.count());

            gcRatioWeightedTotal += cnSegment.median() * sqrtCount;
            fittedDiffWeightedTotalAbs += abs(gcRatioFitDiff) * sqrtCount;
        }

        double weightedResiduals = gcRatioWeightedTotal > 0 ? fittedDiffWeightedTotalAbs / gcRatioWeightedTotal : 0;
        mResiduals = weightedResiduals;

        mEstimatedPurity = 2 * mFitCoefficient / (1 + (2 - samplePloidy) * mFitCoefficient);

        mValid = true;
    }

}
