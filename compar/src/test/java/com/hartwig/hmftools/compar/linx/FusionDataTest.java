package com.hartwig.hmftools.compar.linx;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class FusionDataTest extends ComparableItemTest<FusionData, FusionComparer, TestFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new FusionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestFusionDataBuilder.BUILDER;

        FusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.ReportedType.toString(), b -> b.reportedType = alternateValueSource.Fusion.reportedType());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.Phased.toString(), b -> b.phased = alternateValueSource.Fusion.phased());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.Likelihood.toString(), b -> b.likelihood = alternateValueSource.Fusion.likelihood());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.FusedExonUp.toString(), b -> b.exonUp = alternateValueSource.Fusion.fusedExonUp());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.FusedExonDown.toString(), b -> b.exonDown = alternateValueSource.Fusion.fusedExonDown());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.ChainLinks.toString(), b -> b.chainLinks = alternateValueSource.Fusion.chainLinks());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.ChainTerminated.toString(), b -> b.chainTerminated = alternateValueSource.Fusion.chainTerminated());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.DomainsKept.toString(), b -> b.domainsKept = alternateValueSource.Fusion.domainsKept());
        fieldToAlternateValueInitializer.put(FusionComparer.Fields.DomainsLost.toString(), b -> b.domainsLost = alternateValueSource.Fusion.domainsLost());

        nameToAlternateIndexInitializer = Map.of("FusionName", b -> b.fusionName = alternateValueSource.GeneMappedName);
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FusionComparer.Fields.Reported.toString(), b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
