package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import java.util.List;

import com.hartwig.hmftools.common.sigs.LeastSquaresFit;

public final class CnPurityCalculator
{
    public static CnFitResult calculatePurity(final List<CopyNumberGcData> copyNumberSegments, final double samplePloidy)
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

        double fitCoefficient = lsqFit.getContribs()[0];
        double fitIntercept = avgGcRatio - avgCopyNumber * fitCoefficient;

        double fittedDiffWeightedTotalAbs = 0;
        double gcRatioWeightedTotal = 0;

        for(CopyNumberGcData cnSegment : copyNumberSegments)
        {
            if(!cnSegment.IsValid)
                continue;

            double fittedGcRatio = fitIntercept + cnSegment.CopyNumberLevel * fitCoefficient;
            double gcRatioFitDiff = cnSegment.median() - fittedGcRatio;

            double sqrtCount = sqrt(cnSegment.count());

            gcRatioWeightedTotal += cnSegment.median() * sqrtCount;
            fittedDiffWeightedTotalAbs += abs(gcRatioFitDiff) * sqrtCount;
        }

        double weightedResiduals = gcRatioWeightedTotal > 0 ? fittedDiffWeightedTotalAbs / gcRatioWeightedTotal : 0;
        double residuals = weightedResiduals;

        double estimatedPurity = max(min(2 * fitCoefficient / (1 + (2 - samplePloidy) * fitCoefficient), 1), 0);

        return new CnFitResult(fitCoefficient, fitIntercept, estimatedPurity, residuals);
    }
}
