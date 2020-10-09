package com.hartwig.hmftools.healthchecker;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class HealthCheckEvaluationTest {

    @Test
    public void refCoverageChecksWork() {
        QCValue ref10xSucceed = ImmutableQCValue.of(QCValueType.REF_COVERAGE_10X, "0.95");
        QCValue ref10xFail = ImmutableQCValue.of(QCValueType.REF_COVERAGE_10X, "0.85");

        QCValue ref20xSucceed = ImmutableQCValue.of(QCValueType.REF_COVERAGE_20X, "0.8");
        QCValue ref20xFail = ImmutableQCValue.of(QCValueType.REF_COVERAGE_20X, "0.6");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xSucceed, ref20xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xSucceed, ref20xFail)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xFail, ref20xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xFail, ref20xFail)));
    }

    @Test
    public void tumorCoverageChecksWork() {
        QCValue tumor30xSucceed = ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_30X, "0.9");
        QCValue tumor30xFail = ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_30X, "0.7");

        QCValue tumor60xSucceed = ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_60X, "0.7");
        QCValue tumor60xFail = ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_60X, "0.5");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xSucceed, tumor60xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xSucceed, tumor60xFail)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xFail, tumor60xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xFail, tumor60xFail)));
    }

    @Test
    public void amberMeanBafChecksWork() {
        QCValue amberMeanBafCorrect = ImmutableQCValue.of(QCValueType.AMBER_MEAN_BAF, "0.498");
        QCValue amberMeanBafTooLow = ImmutableQCValue.of(QCValueType.AMBER_MEAN_BAF, "0.4");
        QCValue amberMeanBafTooHigh = ImmutableQCValue.of(QCValueType.AMBER_MEAN_BAF, "0.7");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(amberMeanBafCorrect)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(amberMeanBafTooLow)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(amberMeanBafTooHigh)));
    }

    @Test
    public void amberContaminationCheckWorks() {
        QCValue amberContaminationCorrect = ImmutableQCValue.of(QCValueType.AMBER_CONTAMINATION, "0");
        QCValue amberContaminationCorrectButNonZero = ImmutableQCValue.of(QCValueType.AMBER_CONTAMINATION, "0.01");
        QCValue amberContaminationTooHigh = ImmutableQCValue.of(QCValueType.AMBER_CONTAMINATION, "0.2");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(amberContaminationCorrect)));
        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(amberContaminationCorrectButNonZero)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(amberContaminationTooHigh)));
    }

    @Test
    public void purpleQCCheckWorks() {
        QCValue purpleQCCheckCorrect = ImmutableQCValue.of(QCValueType.PURPLE_QC_STATUS, "PASS");
        QCValue purpleQCCheckFail = ImmutableQCValue.of(QCValueType.PURPLE_QC_STATUS, "FAIL");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleQCCheckCorrect)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleQCCheckFail)));
    }

    @Test
    public void purpleContaminationCheckWorks() {
        QCValue purpleContaminationCorrect = ImmutableQCValue.of(QCValueType.PURPLE_CONTAMINATION, "0");
        QCValue purpleContaminationCorrectButNonZero = ImmutableQCValue.of(QCValueType.PURPLE_CONTAMINATION, "0.01");
        QCValue purpleContaminationTooHigh = ImmutableQCValue.of(QCValueType.PURPLE_CONTAMINATION, "0.2");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationCorrect)));
        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationCorrectButNonZero)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationTooHigh)));
    }
}