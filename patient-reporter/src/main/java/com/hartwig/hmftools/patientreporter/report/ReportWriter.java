package com.hartwig.hmftools.patientreporter.report;

import java.io.FileNotFoundException;

import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public interface ReportWriter {

    @NotNull
    String writeSequenceReport(@NotNull final PatientReport report, @NotNull final Slicer hmfSlicingRegion,
            @NotNull final DrupFilter drupFilter) throws FileNotFoundException, DRException;

    @NotNull
    String writeNonSequenceableReport(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage, @NotNull final NotSequenceableReason reason)
            throws FileNotFoundException, DRException;

}
