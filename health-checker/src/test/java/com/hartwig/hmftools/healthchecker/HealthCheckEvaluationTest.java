package com.hartwig.hmftools.healthchecker;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HealthCheckEvaluationTest {

    @Test
    public void refCoverageChecksWork() {
        QCValue ref10xSucceed = qc(QCValueType.REF_COVERAGE_10X, "0.95");
        QCValue ref10xFail = qc(QCValueType.REF_COVERAGE_10X, "0.85");

        QCValue ref20xSucceed = qc(QCValueType.REF_COVERAGE_20X, "0.8");
        QCValue ref20xFail = qc(QCValueType.REF_COVERAGE_20X, "0.6");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xSucceed, ref20xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xSucceed, ref20xFail)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xFail, ref20xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(ref10xFail, ref20xFail)));
    }

    @Test
    public void tumorCoverageChecksWork() {
        QCValue tumor30xSucceed = qc(QCValueType.TUM_COVERAGE_30X, "0.9");
        QCValue tumor30xFail = qc(QCValueType.TUM_COVERAGE_30X, "0.7");

        QCValue tumor60xSucceed = qc(QCValueType.TUM_COVERAGE_60X, "0.7");
        QCValue tumor60xFail = qc(QCValueType.TUM_COVERAGE_60X, "0.5");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xSucceed, tumor60xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xSucceed, tumor60xFail)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xFail, tumor60xSucceed)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(tumor30xFail, tumor60xFail)));
    }

    @Test
    public void purpleQCCheckWorks() {
        QCValue purpleQCCheckCorrect = qc(QCValueType.PURPLE_QC_STATUS, "PASS");
        QCValue purpleQCCheckFail = qc(QCValueType.PURPLE_QC_STATUS, "FAIL");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleQCCheckCorrect)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleQCCheckFail)));
    }

    @Test
    public void purpleContaminationCheckWorks() {
        QCValue purpleContaminationCorrect = qc(QCValueType.PURPLE_CONTAMINATION, "0");
        QCValue purpleContaminationCorrectButNonZero = qc(QCValueType.PURPLE_CONTAMINATION, "0.01");
        QCValue purpleContaminationTooHigh = qc(QCValueType.PURPLE_CONTAMINATION, "0.2");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationCorrect)));
        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationCorrectButNonZero)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(purpleContaminationTooHigh)));
    }

    @Test
    public void refFlagstatProportionMappedCheckWorks() {
        QCValue proportionMappedCorrect = qc(QCValueType.REF_PROPORTION_MAPPED, "0.97");
        QCValue proportionMappedTooLow = qc(QCValueType.REF_PROPORTION_MAPPED, "0.4");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(proportionMappedCorrect)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(proportionMappedTooLow)));
    }

    @Test
    public void tumFlagstatProportionMappedCheckWorks() {
        QCValue proportionMappedCorrect = qc(QCValueType.TUM_PROPORTION_MAPPED, "0.97");
        QCValue proportionMappedTooLow = qc(QCValueType.TUM_PROPORTION_MAPPED, "0.4");

        assertTrue(HealthCheckEvaluation.isPass(Lists.newArrayList(proportionMappedCorrect)));
        assertFalse(HealthCheckEvaluation.isPass(Lists.newArrayList(proportionMappedTooLow)));
    }

    @NotNull
    private static QCValue qc(@NotNull QCValueType type, @NotNull String value) {
        return ImmutableQCValue.builder().type(type).value(value).build();
    }
}