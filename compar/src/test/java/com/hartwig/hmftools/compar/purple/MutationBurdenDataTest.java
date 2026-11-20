package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.HashMap;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class MutationBurdenDataTest extends ComparableItemTest<MutationBurdenData, MutationBurdenComparer, TestPurityDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new MutationBurdenComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestPurityDataBuilder.MUTATION_BURDEN_BUILDER;
        PurityContext alternateValues = builder.createWithAlternateDefaults().Purity;

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.TmbPerMb.toString(), b -> b.tmb = alternateValues.tumorMutationalBurdenPerMb());
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.MsIndelsPerMb.toString(), b -> b.msIndels = alternateValues.microsatelliteIndelsPerMb());
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.Tml.toString(), b -> b.tml = alternateValues.tumorMutationalLoad());
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.MsStatus.toString(), b -> b.msStatus = alternateValues.microsatelliteStatus());
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.TmbStatus.toString(), b -> b.tmbStatus = alternateValues.tumorMutationalBurdenStatus());
        fieldToAlternateValueInitializer.put(
                MutationBurdenComparer.Fields.TmlStatus.toString(), b -> b.tmlStatus = alternateValues.tumorMutationalLoadStatus());

        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
