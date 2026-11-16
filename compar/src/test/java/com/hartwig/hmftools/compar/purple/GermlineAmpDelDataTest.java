package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class GermlineAmpDelDataTest
        extends ComparableItemTest<GermlineAmpDelData, GermlineAmpDelComparer, TestGermlineAmpDelDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GermlineAmpDelComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestGermlineAmpDelDataBuilder.BUILDER;
        GermlineAmpDelData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                GermlineAmpDelComparer.Fields.GermlineStatus.toString(), b -> b.germlineStatus = alternateValueSource.AmpDelData.NormalStatus,
                GermlineAmpDelComparer.Fields.TumorStatus.toString(), b -> b.tumorStatus = alternateValueSource.AmpDelData.TumorStatus,
                GermlineAmpDelComparer.Fields.GermlineCopyNumber.toString(), b -> b.germlineCopyNumber = alternateValueSource.AmpDelData.GermlineCopyNumber,
                GermlineAmpDelComparer.Fields.TumorCopyNumber.toString(), b -> b.tumorCopyNumber = alternateValueSource.AmpDelData.TumorCopyNumber
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.AmpDelData.GeneName);
        reportabilityFieldToFalseReportabilityInitializer = Map.of(GermlineAmpDelComparer.Fields.Reported.toString(), b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
