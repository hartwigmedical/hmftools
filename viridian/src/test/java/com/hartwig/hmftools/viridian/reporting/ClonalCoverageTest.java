package com.hartwig.hmftools.viridian.reporting;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.OptionalDouble;
import java.util.Set;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.BamMetricsTestFactory;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricSummary;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;

import org.junit.Test;

public class ClonalCoverageTest
{
    private static final double EPSILON = 1e-9;

    @Test
    public void testSimpleFullPurityCase()
    {
        // rho=1, psi=2: denom = 1*2 + 2*0 = 2; depth per copy = D_host*rho/denom = 100*1/2 = 50
        OptionalDouble depthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(
                purityContext(1.0, 2.0, FittedPurityMethod.SOMATIC, Set.of(PurpleQCStatus.PASS)),
                bamMetrics(100.0));

        assertTrue(depthPerCopy.isPresent());
        assertEquals(50.0, depthPerCopy.getAsDouble(), EPSILON);
    }

    @Test
    public void testNormalisesForSubclonalPurity()
    {
        // denom = rho*psi + 2(1-rho) = 0.5*2 + 2*0.5 = 2; depth per copy = 30*0.5/2 = 7.5
        OptionalDouble depthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(
                purityContext(0.5, 2.0, FittedPurityMethod.SOMATIC, Set.of(PurpleQCStatus.PASS)),
                bamMetrics(30.0));

        assertTrue(depthPerCopy.isPresent());
        assertEquals(7.5, depthPerCopy.getAsDouble(), EPSILON);
    }

    @Test
    public void testComputesWhenFitOnlyWarns()
    {
        OptionalDouble depthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(
                purityContext(0.8, 2.0, FittedPurityMethod.SOMATIC, Set.of(PurpleQCStatus.WARN_DELETED_GENES)),
                bamMetrics(30.0));

        assertTrue(depthPerCopy.isPresent());
    }

    @Test
    public void testEmptyWhenPurityFitFailed()
    {
        OptionalDouble depthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(
                purityContext(0.5, 2.0, FittedPurityMethod.SOMATIC, Set.of(PurpleQCStatus.FAIL_CONTAMINATION)),
                bamMetrics(30.0));

        assertFalse(depthPerCopy.isPresent());
    }

    @Test
    public void testEmptyForNoTumorFit()
    {
        OptionalDouble depthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(
                purityContext(0.0, 2.0, FittedPurityMethod.NO_TUMOR, Set.of(PurpleQCStatus.FAIL_NO_TUMOR)),
                bamMetrics(30.0));

        assertFalse(depthPerCopy.isPresent());
    }

    private static BamMetricSummary bamMetrics(double meanCoverage)
    {
        return ImmutableBamMetricSummary.builder()
                .from(BamMetricsTestFactory.createMinimalTestWGSMetrics())
                .meanCoverage(meanCoverage)
                .build();
    }

    private static PurityContext purityContext(
            double purity, double ploidy, FittedPurityMethod method, Set<PurpleQCStatus> qcStatus)
    {
        FittedPurity bestFit = new FittedPurity(purity, -1, ploidy, -1, -1, -1);
        PurpleQC qc = ImmutablePurpleQC.builder()
                .status(qcStatus)
                .method(method)
                .copyNumberSegments(-1)
                .unsupportedCopyNumberSegments(-1)
                .purity(purity)
                .contamination(0)
                .cobaltGender(Gender.MALE)
                .amberGender(Gender.MALE)
                .deletedGenes(-1)
                .amberMeanDepth(-1)
                .lohPercent(-1)
                .tincLevel(0)
                .build();
        return ImmutablePurityContext.builder()
                .gender(Gender.MALE)
                .bestFit(bestFit)
                .method(method)
                .qc(qc)
                .microsatelliteIndelsPerMb(0)
                .tumorMutationalBurdenPerMb(0)
                .tumorMutationalLoad(0)
                .svTumorMutationalBurden(0)
                .microsatelliteStatus(MicrosatelliteStatus.MSS)
                .tumorMutationalLoadStatus(TumorMutationalStatus.LOW)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.LOW)
                .runMode(RunMode.TUMOR_GERMLINE)
                .targeted(false)
                .score(ImmutableFittedPurityScore.builder()
                        .minPurity(-1).maxPurity(-1).minPloidy(-1).maxPloidy(-1)
                        .minDiploidProportion(-1).maxDiploidProportion(-1)
                        .build())
                .polyClonalProportion(-1)
                .wholeGenomeDuplication(false)
                .build();
    }
}
