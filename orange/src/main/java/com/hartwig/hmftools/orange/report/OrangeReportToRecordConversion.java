package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.orange.algo.OrangeReport;

public class OrangeReportToRecordConversion {

    private OrangeReportToRecordConversion() {
    }

    public static OrangeRecord convert(OrangeReport report) {
        return ImmutableOrangeRecord.builder()
                .sampleId(report.sampleId())
                .experimentDate(report.experimentDate())
                .refGenomeVersion(report.refGenomeVersion())
                .purple(report.purple())
                .linx(report.linx())
                .lilac(report.lilac())
                .virusInterpreter(report.virusInterpreter())
                .chord(report.chord())
                .cuppa(report.cuppa())
                .peach(report.peach())
                .plots(report.plots())
                .build();
    }
}
