package com.hartwig.hmftools.orange.algo.viridian;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.viridian.ClonalityStatus;
import com.hartwig.hmftools.datamodel.viridian.ImmutableViridianRecord;
import com.hartwig.hmftools.datamodel.viridian.ImmutableVirusReport;
import com.hartwig.hmftools.datamodel.viridian.ReportQcStatus;
import com.hartwig.hmftools.datamodel.viridian.ViridianRecord;
import com.hartwig.hmftools.datamodel.viridian.VirusInterpretation;
import com.hartwig.hmftools.datamodel.viridian.VirusReport;

public final class ViridianInterpreter
{
    public static ViridianRecord interpret(final List<com.hartwig.hmftools.common.viridian.VirusReport> reports)
    {
        List<VirusReport> viruses = reports.stream()
                .map(report -> ImmutableVirusReport.builder()
                        .name(report.virusName())
                        .virusGroup(report.virusGroup())
                        .qcStatus(ReportQcStatus.valueOf(report.qcStatus().name()))
                        .integrations(report.integrations())
                        .interpretation(report.reportingType() != null
                                ? VirusInterpretation.valueOf(report.reportingType().name())
                                : null)
                        .percentageCovered(report.coveragePercentage())
                        .meanCoverage(report.meanDepth())
                        .reported(report.isReported())
                        .clonalityStatus(ClonalityStatus.valueOf(report.clonalityStatus().name()))
                        .copiesPerTumorCell(report.copiesPerTumorCell())
                        .build())
                .collect(Collectors.toList());

        return ImmutableViridianRecord.builder()
                .allViruses(viruses)
                .reportableViruses(viruses.stream().filter(VirusReport::reported).toList())
                .build();
    }
}
