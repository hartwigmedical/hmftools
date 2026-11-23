package com.hartwig.hmftools.compar.virus;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.viridian.ClonalityStatus;
import com.hartwig.hmftools.common.viridian.ReportQcStatus;
import com.hartwig.hmftools.common.viridian.VirusReport;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestVirusDataBuilder
{
    public String virusName = "Human papillomavirus type 16";
    public String virusGroup = "Human papillomavirus";
    public VirusType reportingType = VirusType.HPV;
    public int integrations = 4;
    public double coverage = 0.9;
    public double meanDepth = 20;
    public Double copiesPerTumorCell = 2.0;
    public ClonalityStatus clonalityStatus = ClonalityStatus.CLONAL;
    public boolean reported = true;

    private static final Consumer<TestVirusDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.virusName = "Human gammaherpesvirus 4";
        b.virusGroup = "Human gammaherpesvirus 4";
        b.reportingType = VirusType.EBV;
        b.integrations = 40;
        b.coverage = 0.2;
        b.meanDepth = 2;
        b.copiesPerTumorCell = 0.2;
        b.clonalityStatus = ClonalityStatus.NOT_CLONAL;
        b.reported = false;
    };

    public static final TestComparableItemBuilder<TestVirusDataBuilder, VirusData> BUILDER =
            new TestComparableItemBuilder<>(TestVirusDataBuilder::new, TestVirusDataBuilder::build, ALTERNATE_INITIALIZER);

    private VirusData build()
    {
        return VirusData.fromViridian(
                new VirusReport(
                        virusName, virusGroup, ReportQcStatus.NO_ABNORMALITIES, reportingType, integrations, coverage, meanDepth,
                        copiesPerTumorCell, clonalityStatus, reported),
                new VirusComparer(null, Collections.emptyMap()).fieldsList());
    }
}
