package com.hartwig.hmftools.compar.peach;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class PeachDataTest extends ComparableItemTest<PeachData, PeachComparer, TestPeachDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new PeachComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestPeachDataBuilder.BUILDER;
        PeachData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                PeachComparer.Fields.AlleleCount.toString(), b -> b.alleleCount = alternateValueSource.Genotype.alleleCount(),
                PeachComparer.Fields.Function.toString(), b -> b.function = alternateValueSource.Genotype.function(),
                PeachComparer.Fields.Drugs.toString(), b -> b.linkedDrugs = alternateValueSource.Genotype.linkedDrugs(),
                PeachComparer.Fields.PrescriptionUrls.toString(), b -> b.prescriptionUrls = alternateValueSource.Genotype.urlPrescriptionInfo()
        );
        nameToAlternateIndexInitializer = Map.of(
                "Gene", b -> b.gene = alternateValueSource.Genotype.gene(),
                "Allele", b -> b.allele = alternateValueSource.Genotype.allele()
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
