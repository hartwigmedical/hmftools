package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_DRIVER_LIKELIHOOD;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_INTEGRATIONS;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_MEAN_COVERAGE;

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
        comparer = new VirusComparer(new ComparConfig());
        builder = TestVirusDataBuilder.BUILDER;
        VirusData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_INTEGRATIONS, b -> b.integrations = alternateValueSource.Virus.integrations(),
                FLD_MEAN_COVERAGE, b -> b.meanCoverage = alternateValueSource.Virus.meanCoverage(),
                FLD_DRIVER_LIKELIHOOD, b -> b.driverLikelihood = alternateValueSource.Virus.virusDriverLikelihoodType()
        );
        nameToAlternateIndexInitializer = Map.of("name", b -> b.name = alternateValueSource.Virus.name());
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
