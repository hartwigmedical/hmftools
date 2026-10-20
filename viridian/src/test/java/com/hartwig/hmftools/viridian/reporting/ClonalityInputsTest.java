package com.hartwig.hmftools.viridian.reporting;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

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

public class ClonalityInputsTest
{
    @Test
    public void testSimpleFullPurityCase()
    {
        // rho=1, psi=2: denom = 1*2 + 2*0 = 2; depth per copy = D_host*rho/denom = 100*1/2 = 50
        assertEquals(
                new ClonalityInputs(50.0, null),
                clonalityInputs(1.0, FittedPurityMethod.SOMATIC, 100.0, PurpleQCStatus.PASS));
    }

    @Test
    public void testNormalisesForSubclonalPurity()
    {
        // denom = rho*psi + 2(1-rho) = 0.5*2 + 2*0.5 = 2; depth per copy = 30*0.5/2 = 7.5
        assertEquals(
                new ClonalityInputs(7.5, null),
                clonalityInputs(0.5, FittedPurityMethod.SOMATIC, 30.0, PurpleQCStatus.PASS));
    }

    @Test
    public void testComputesWhenFitOnlyWarns()
    {
        assertEquals(
                new ClonalityInputs(15.0, null),
                clonalityInputs(1.0, FittedPurityMethod.SOMATIC, 30.0, PurpleQCStatus.WARN_DELETED_GENES));
    }

    // Only no tumor and contamination block the calculation; other Purple failures are allowed.
    @Test
    public void testComputesWhenTumorInNormalFails()
    {
        assertEquals(
                new ClonalityInputs(15.0, null),
                clonalityInputs(1.0, FittedPurityMethod.SOMATIC, 30.0, PurpleQCStatus.FAIL_TINC));
    }

    @Test
    public void testContaminationFailure()
    {
        assertEquals(
                new ClonalityInputs(null, PurpleQCStatus.FAIL_CONTAMINATION),
                clonalityInputs(0.5, FittedPurityMethod.SOMATIC, 30.0, PurpleQCStatus.FAIL_CONTAMINATION));
    }

    @Test
    public void testNoTumorFailure()
    {
        assertEquals(
                new ClonalityInputs(null, PurpleQCStatus.FAIL_NO_TUMOR),
                clonalityInputs(0.08, FittedPurityMethod.NO_TUMOR, 30.0, PurpleQCStatus.FAIL_NO_TUMOR));
    }

    // With no tumor, copies per tumor cell is meaningless whatever the contamination.
    @Test
    public void testNoTumorWinsOverContamination()
    {
        assertEquals(
                new ClonalityInputs(null, PurpleQCStatus.FAIL_NO_TUMOR),
                clonalityInputs(
                        0.08, FittedPurityMethod.NO_TUMOR, 30.0,
                        PurpleQCStatus.FAIL_NO_TUMOR, PurpleQCStatus.FAIL_CONTAMINATION));
    }

    // Exactly one of the depth per copy and a clonality-blocking Purple failure.
    @Test
    public void testRejectsInvalidCombinations()
    {
        assertThrows(IllegalArgumentException.class, () -> new ClonalityInputs(null, null));
        assertThrows(IllegalArgumentException.class, () -> new ClonalityInputs(15.0, PurpleQCStatus.FAIL_NO_TUMOR));
        assertThrows(IllegalArgumentException.class, () -> new ClonalityInputs(null, PurpleQCStatus.FAIL_TINC));
    }

    private static ClonalityInputs clonalityInputs(
            double purity, FittedPurityMethod method, double hostMeanDepth, PurpleQCStatus... qcStatus)
    {
        return ClonalityInputs.from(purityContext(purity, 2.0, method, Set.of(qcStatus)), bamMetrics(hostMeanDepth));
    }

    private static BamMetricSummary bamMetrics(double meanDepth)
    {
        return ImmutableBamMetricSummary.builder()
                .from(BamMetricsTestFactory.createMinimalTestWGSMetrics())
                .meanCoverage(meanDepth)
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
