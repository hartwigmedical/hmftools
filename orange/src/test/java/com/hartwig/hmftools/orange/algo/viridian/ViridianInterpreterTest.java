package com.hartwig.hmftools.orange.algo.viridian;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.viridian.ClonalityStatus;
import com.hartwig.hmftools.common.viridian.ReportQcStatus;
import com.hartwig.hmftools.common.viridian.VirusReport;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.datamodel.viridian.ImmutableViridianRecord;
import com.hartwig.hmftools.datamodel.viridian.ImmutableVirusReport;
import com.hartwig.hmftools.datamodel.viridian.VirusInterpretation;

import org.junit.Test;

public class ViridianInterpreterTest
{
    @Test
    public void convertsEachVirusConstantProperly()
    {
        for(final VirusType value : VirusType.values())
        {
            VirusInterpretation.valueOf(value.name());
        }
        for(final ReportQcStatus value : ReportQcStatus.values())
        {
            com.hartwig.hmftools.datamodel.viridian.ReportQcStatus.valueOf(value.name());
        }
        for(final ClonalityStatus value : ClonalityStatus.values())
        {
            com.hartwig.hmftools.datamodel.viridian.ClonalityStatus.valueOf(value.name());
        }
    }

    @Test
    public void testInterpret()
    {
        VirusReport reported = new VirusReport(
                "Human papillomavirus type 16", "Human papillomavirus", ReportQcStatus.AMBIGUOUS_VIRUS, VirusType.HPV,
                2, 90.0, 20.0, 2.5, ClonalityStatus.CLONAL, true);
        VirusReport notReported = new VirusReport(
                "JC polyomavirus", "JC polyomavirus", ReportQcStatus.NO_ABNORMALITIES, null,
                0, 40.0, 1.0, null, ClonalityStatus.PURPLE_FAIL_NO_TUMOR, false);

        com.hartwig.hmftools.datamodel.viridian.VirusReport expectedReported = ImmutableVirusReport.builder()
                .name("Human papillomavirus type 16")
                .virusGroup("Human papillomavirus")
                .qcStatus(com.hartwig.hmftools.datamodel.viridian.ReportQcStatus.AMBIGUOUS_VIRUS)
                .integrations(2)
                .interpretation(VirusInterpretation.HPV)
                .percentageCovered(90.0)
                .meanCoverage(20.0)
                .reported(true)
                .clonalityStatus(com.hartwig.hmftools.datamodel.viridian.ClonalityStatus.CLONAL)
                .copiesPerTumorCell(2.5)
                .build();
        com.hartwig.hmftools.datamodel.viridian.VirusReport expectedNotReported = ImmutableVirusReport.builder()
                .name("JC polyomavirus")
                .virusGroup("JC polyomavirus")
                .qcStatus(com.hartwig.hmftools.datamodel.viridian.ReportQcStatus.NO_ABNORMALITIES)
                .integrations(0)
                .interpretation(null)
                .percentageCovered(40.0)
                .meanCoverage(1.0)
                .reported(false)
                .clonalityStatus(com.hartwig.hmftools.datamodel.viridian.ClonalityStatus.PURPLE_FAIL_NO_TUMOR)
                .copiesPerTumorCell(null)
                .build();

        assertEquals(
                ImmutableViridianRecord.builder()
                        .allViruses(List.of(expectedReported, expectedNotReported))
                        .reportableViruses(List.of(expectedReported))
                        .build(),
                ViridianInterpreter.interpret(List.of(reported, notReported)));
    }
}
