package com.hartwig.hmftools.compar.virus;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class VirusDataTest extends ComparableItemTest<VirusData, VirusComparer, TestVirusDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new VirusComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestVirusDataBuilder.BUILDER;
        VirusData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                VirusComparer.Fields.Integrations.toString(), b -> b.integrations = alternateValueSource.Virus.integrations(),
                VirusComparer.Fields.MeanCoverage.toString(), b -> b.meanCoverage = alternateValueSource.Virus.meanCoverage(),
                VirusComparer.Fields.DriverLikelihood.toString(), b -> b.driverLikelihood = alternateValueSource.Virus.virusDriverLikelihoodType()
        );
        nameToAlternateIndexInitializer = Map.of("name", b -> b.name = alternateValueSource.Virus.name());
        reportabilityFieldToFalseReportabilityInitializer = Map.of(VirusComparer.Fields.Reported.toString(), b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
