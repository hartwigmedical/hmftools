package com.hartwig.hmftools.patientreporter.util;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

public final class FindingsToCSV {

    private FindingsToCSV() {
    }

    @NotNull
    public static List<String> varToCSV(@NotNull final List<VariantReport> reports) {
        final List<String> lines = Lists.newArrayList();
        lines.add("GENE,POSITION,REF,ALT,TRANSCRIPT,CDS,P,CONSEQUENCE,COSMIC_ID,ALLELE_READ_COUNT,TOTAL_READ_COUNT");
        lines.addAll(reports.stream().map(
                report -> report.gene() + "," + report.position() + "," + report.ref() + "," + report.alt() + ","
                        + report.transcript() + "," + report.hgvsCoding() + "," + report.hgvsProtein() + ","
                        + report.consequence() + "," + report.cosmicID() + "," + Integer.toString(
                        report.alleleReadCount()) + "," + Integer.toString(report.totalReadCount())).
                collect(Collectors.toList()));
        return lines;
    }

    @NotNull
    public static List<String> cnvToCSV(@NotNull final List<CopyNumberReport> reports) {
        final List<String> lines = Lists.newArrayList();
        lines.add("GENE,TRANSCRIPT,FINDING");
        lines.addAll(reports.stream().map(report -> report.gene() + "," + report.transcript() + "," + Integer.toString(
                report.copyNumber())).collect(Collectors.toList()));
        return lines;
    }
}
