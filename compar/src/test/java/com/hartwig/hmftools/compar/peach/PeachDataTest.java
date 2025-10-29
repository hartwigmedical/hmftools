package com.hartwig.hmftools.compar.peach;

import static com.hartwig.hmftools.compar.peach.PeachData.FLD_ALLELE_COUNT;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_DRUGS;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_FUNCTION;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_PRESCRIPTION_URLS;

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
        comparer = new PeachComparer(new ComparConfig());
        builder = TestPeachDataBuilder.BUILDER;
        PeachData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_ALLELE_COUNT, b -> b.alleleCount = alternateValueSource.Genotype.alleleCount(),
                FLD_FUNCTION, b -> b.function = alternateValueSource.Genotype.function(),
                FLD_DRUGS, b -> b.linkedDrugs = alternateValueSource.Genotype.linkedDrugs(),
                FLD_PRESCRIPTION_URLS, b -> b.prescriptionUrls = alternateValueSource.Genotype.urlPrescriptionInfo()
        );
        nameToAlternateIndexInitializer = Map.of(
                "Gene", b -> b.gene = alternateValueSource.Genotype.gene(),
                "Allele", b -> b.allele = alternateValueSource.Genotype.allele()
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
