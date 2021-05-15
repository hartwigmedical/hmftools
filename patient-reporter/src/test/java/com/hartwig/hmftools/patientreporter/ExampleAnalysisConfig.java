package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ExampleAnalysisConfig {

    @NotNull
    private final String sampleId;
    private final boolean isCorrectionReport;
    @Nullable
    private final String comments;
    @NotNull
    private final QsFormNumber qcForNumber;
    private final boolean hasReliablePurity;
    private final double impliedTumorPurity;
    private final boolean includeSummary;
    private final boolean reportGermline;
    @NotNull
    private final LimsCohortConfig limsCohortConfig;

    private ExampleAnalysisConfig(@NotNull final String sampleId, final boolean isCorrectionReport, @Nullable final String comments,
            @NotNull final QsFormNumber qcForNumber, final boolean hasReliablePurity, final double impliedTumorPurity,
            final boolean includeSummary, final boolean reportGermline, @NotNull final LimsCohortConfig limsCohortConfig) {
        this.sampleId = sampleId;
        this.isCorrectionReport = isCorrectionReport;
        this.comments = comments;
        this.qcForNumber = qcForNumber;
        this.hasReliablePurity = hasReliablePurity;
        this.impliedTumorPurity = impliedTumorPurity;
        this.includeSummary = includeSummary;
        this.reportGermline = reportGermline;
        this.limsCohortConfig = limsCohortConfig;
    }

    @NotNull
    public String sampleId() {
        return sampleId;
    }

    public boolean isCorrectionReport() {
        return isCorrectionReport;
    }

    @Nullable
    public String comments() {
        return comments;
    }

    @NotNull
    public QsFormNumber qcForNumber() {
        return qcForNumber;
    }

    public boolean hasReliablePurity() {
        return hasReliablePurity;
    }

    public double impliedTumorPurity() {
        return impliedTumorPurity;
    }

    public boolean includeSummary() {
        return includeSummary;
    }

    public boolean reportGermline() {
        return reportGermline;
    }

    @NotNull
    public LimsCohortConfig limsCohortConfig() {
        return limsCohortConfig;
    }

    public static class Builder {
        @NotNull
        private String sampleId = "COLO829T";
        private boolean isCorrectionReport = false;
        @Nullable
        private String comments = null;
        @NotNull
        private QsFormNumber qcForNumber = QsFormNumber.FOR_080;
        private boolean hasReliablePurity = true;
        private double impliedTumorPurity = 1D;
        private boolean includeSummary = true;
        private boolean reportGermline = false;
        @NotNull
        private LimsCohortConfig limsCohortConfig = LimsCohortTestFactory.createCOLOCohortConfig();

        public Builder() {
        }

        @NotNull
        public Builder sampleId(@NotNull final String sampleId) {
            this.sampleId = sampleId;
            return this;
        }

        @NotNull
        public Builder isCorrectionReport(final boolean isCorrectionReport) {
            this.isCorrectionReport = isCorrectionReport;
            return this;
        }

        @NotNull
        public Builder comments(@Nullable final String comments) {
            this.comments = comments;
            return this;
        }

        @NotNull
        public Builder qcForNumber(@NotNull final QsFormNumber qcForNumber) {
            this.qcForNumber = qcForNumber;
            return this;
        }

        @NotNull
        public Builder hasReliablePurity(final boolean hasReliablePurity) {
            this.hasReliablePurity = hasReliablePurity;
            return this;
        }

        @NotNull
        public Builder impliedTumorPurity(final double impliedTumorPurity) {
            this.impliedTumorPurity = impliedTumorPurity;
            return this;
        }

        @NotNull
        public Builder includeSummary(final boolean includeSummary) {
            this.includeSummary = includeSummary;
            return this;
        }

        @NotNull
        public Builder reportGermline(final boolean reportGermline) {
            this.reportGermline = reportGermline;
            return this;
        }

        @NotNull
        public Builder limsCohortConfig(@NotNull final LimsCohortConfig limsCohortConfig) {
            this.limsCohortConfig = limsCohortConfig;
            return this;
        }

        @NotNull
        public ExampleAnalysisConfig build() {
            return new ExampleAnalysisConfig(sampleId, isCorrectionReport,
                    comments,
                    qcForNumber,
                    hasReliablePurity,
                    impliedTumorPurity,
                    includeSummary,
                    reportGermline,
                    limsCohortConfig);
        }
    }
}
