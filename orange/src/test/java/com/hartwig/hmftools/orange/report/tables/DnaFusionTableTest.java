package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.algo.linx.TestLinxRecordFactory.linxRecordFusionBuilder;

import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.hartwig.hmftools.datamodel.common.ImmutableAllelicDepth;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.junit.Test;

public class DnaFusionTableTest
{
    @Test
    public void canBuildTableWithMixedRnaSupport()
    {
        LinxFusion fusionWithRnaSupport = linxRecordFusionBuilder()
                .rnaSupport(ImmutableAllelicDepth.builder().alleleReadCount(13).totalReadCount(20).build())
                .build();
        LinxFusion fusionWithoutRnaSupport = linxRecordFusionBuilder()
                .geneDown("GENE03")
                .rnaSupport(null)
                .build();

        assertNotNull(DnaFusionTable.build(
                "DNA fusions", 100, List.of(fusionWithRnaSupport, fusionWithoutRnaSupport), ReportResources.create()));
    }
}
