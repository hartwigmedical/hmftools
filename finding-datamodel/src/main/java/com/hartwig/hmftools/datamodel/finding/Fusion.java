package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record Fusion(
        @NotNull DriverFields driver,
        @NotNull String geneStart,
        @NotNull String geneContextStart,
        @NotNull String geneTranscriptStart,
        @NotNull String geneEnd,
        @NotNull String geneContextEnd,
        @NotNull String geneTranscriptEnd,
        @NotNull LinxFusionType reportedType,
        @NotNull List<LinxUnreportableReason> unreportedReasons,
        @NotNull FusionPhasedType phased,
        int fusedExonUp,
        int fusedExonDown,
        int chainLinks,
        boolean chainTerminated,
        @NotNull String domainsKept,
        @NotNull String domainsLost,
        double junctionCopyNumber
) implements Driver {

    @NotNull @Override public String findingKey() { return driver.findingKey(); }
    @NotNull @Override public DriverSource driverSource() { return driver.driverSource(); }
    @NotNull @Override public ReportedStatus reportedStatus() { return driver.reportedStatus(); }
    @NotNull @Override public DriverInterpretation driverInterpretation() { return driver.driverInterpretation(); }

    @NotNull
    public String display() {
        return String.format("%s::%s", geneStart, geneEnd);
    }
}
