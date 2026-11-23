package com.hartwig.hmftools.qsee.common;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;

public class ReportableGenes
{
    public static List<String> forSampleType(List<DriverGene> driverGenes, SampleType sampleType)
    {
        List<String> reportableGenes = new ArrayList<>();

        for(DriverGene driverGene : driverGenes)
        {
            boolean isRelevantReportableGene = sampleType == SampleType.TUMOR
                    ? driverGene.reportSomatic()
                    : driverGene.reportGermline();

            if(isRelevantReportableGene)
            {
                reportableGenes.add(driverGene.gene());
            }
        }

        return reportableGenes;
    }
}
