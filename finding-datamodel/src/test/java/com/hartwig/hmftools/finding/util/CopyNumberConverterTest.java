package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import org.junit.Test;

public class CopyNumberConverterTest  {

    @Test
    public void canFilterDisruptions() {
        Disruption disruption = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberConverter.convertDisruption(List.of(disruption)).get(0).reportedStatus());

        Disruption disruptionFilter = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberConverter.convertDisruption(List.of(disruptionFilter)).get(0).reportedStatus());
    }

    @Test
    public void canFilterFusion() {
        Fusion fusion = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberConverter.convertFusion(List.of(fusion)).get(0).reportedStatus());

        Fusion fusionFilter = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberConverter.convertFusion(List.of(fusionFilter)).get(0).reportedStatus());
    }

}