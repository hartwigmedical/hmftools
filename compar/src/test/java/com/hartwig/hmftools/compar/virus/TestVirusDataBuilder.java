package com.hartwig.hmftools.compar.virus;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestVirusDataBuilder
{
    public String name = "Human betaherpesvirus 6B";
    public boolean reported = true;
    public int integrations = 4;
    public double meanCoverage = 20;
    public VirusLikelihoodType driverLikelihood = VirusLikelihoodType.HIGH;

    private static final Consumer<TestVirusDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.name = "ref";
        b.reported = false;
        b.integrations = 40;
        b.meanCoverage = 2;
        b.driverLikelihood = VirusLikelihoodType.LOW;
    };

    public static final TestComparableItemBuilder<TestVirusDataBuilder, VirusData> BUILDER =
            new TestComparableItemBuilder<>(TestVirusDataBuilder::new, TestVirusDataBuilder::build, ALTERNATE_INITIALIZER);

    private VirusData build()
    {
        return new VirusData(
                ImmutableAnnotatedVirus.builder()
                        .name(name)
                        .integrations(integrations)
                        .meanCoverage(meanCoverage)
                        .reported(reported)
                        .virusDriverLikelihoodType(driverLikelihood)
                        .taxid(-1)
                        .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                        .interpretation(VirusType.EBV)
                        .percentageCovered(-1)
                        .expectedClonalCoverage(null)
                        .blacklisted(null)
                        .build()
        );
    }
}
