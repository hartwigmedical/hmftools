package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

import org.apache.logging.log4j.util.Strings;

public class TestGeneCopyNumberDataBuilder
{
    public String gene = "BRAF";
    public double minCopyNumber = 1;
    public double maxCopyNumber = 2;
    public DriverType driverType = DriverType.AMP;
    public ReportedStatus reportedStatus = ReportedStatus.REPORTED;
    public String transcriptId = Strings.EMPTY;
    public boolean isCanonical = true;

    private static final Consumer<TestGeneCopyNumberDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.gene = "BRCA2";
        b.minCopyNumber = 0;
        b.maxCopyNumber = 5;
        b.driverType = DriverType.DEL;
    };

    public static final TestComparableItemBuilder<TestGeneCopyNumberDataBuilder, GeneCopyNumberData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGeneCopyNumberDataBuilder::new,
                    TestGeneCopyNumberDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GeneCopyNumberData build()
    {
        GeneCopyNumber geneCopyNumber = new GeneCopyNumber("", 0, 0, gene, transcriptId, isCanonical,
                    Strings.EMPTY, maxCopyNumber , minCopyNumber, 0, 1, 1,
                    0, 0, 0, 1.0,
                    SegmentSupport.BND, SegmentSupport.BND, CopyNumberMethod.UNKNOWN, minCopyNumber / 2);

        geneCopyNumber.setDriverType(driverType);
        geneCopyNumber.setReportedStatus(reportedStatus);

        return new GeneCopyNumberData(geneCopyNumber, new GeneCopyNumberComparer(null, Collections.emptyMap()).fieldsList());
    }
}
