package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;

public class OrangeReportToRecordConversion {

    private OrangeReportToRecordConversion() {
    }

    public static OrangeRecord convert(OrangeReport report) {
        RefGenomeVersion refGenomeVersion = report.refGenomeVersion();
        return ImmutableOrangeRecord.builder()
                .sampleId(report.sampleId())
                .experimentDate(report.experimentDate())
                .refGenomeVersion(OrangeRefGenomeVersion.valueOf(refGenomeVersion.name()))
                .purple(convert(report.purple()))
                .linx(report.linx())
                .lilac(report.lilac())
                .virusInterpreter(report.virusInterpreter())
                .chord(report.chord())
                .cuppa(report.cuppa())
                .peach(report.peach())
                .plots(report.plots())
                .build();
    }

    private static PurpleRecord convert(PurpleInterpretedData purpleInterpretedData) {
        return ImmutablePurpleRecord.builder()
                .fit(purpleInterpretedData.fit())
                .characteristics(purpleInterpretedData.characteristics())
                .somaticDrivers(purpleInterpretedData.somaticDrivers())
                .germlineDrivers(purpleInterpretedData.germlineDrivers())
                .allSomaticVariants(purpleInterpretedData.allSomaticVariants())
                .reportableSomaticVariants(purpleInterpretedData.reportableSomaticVariants())
                .allGermlineVariants(purpleInterpretedData.allGermlineVariants())
                .reportableGermlineVariants(purpleInterpretedData.reportableGermlineVariants())
                .allSomaticCopyNumbers(purpleInterpretedData.allSomaticCopyNumbers())
                .allSomaticGeneCopyNumbers(purpleInterpretedData.allSomaticGeneCopyNumbers())
                .suspectGeneCopyNumbersWithLOH(purpleInterpretedData.suspectGeneCopyNumbersWithLOH())
                .allSomaticGainsLosses(purpleInterpretedData.reportableSomaticGainsLosses())
                .reportableSomaticGainsLosses(purpleInterpretedData.reportableSomaticGainsLosses())
                .build();
    }
}
