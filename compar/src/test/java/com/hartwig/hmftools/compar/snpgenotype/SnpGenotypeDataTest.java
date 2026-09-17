package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeComparer.FLD_GENOTYPE;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeComparer.FLD_VCF_SAMPLE_ID;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class SnpGenotypeDataTest extends ComparableItemTest<SnpGenotypeData, SnpGenotypeComparer, TestSnpGenotypeDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new SnpGenotypeComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestSnpGenotypeDataBuilder.BUILDER;
        SnpGenotypeData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_ALT, b -> b.alt = alternateValueSource.Alt,
                FLD_GENOTYPE, b -> b.genotype = alternateValueSource.Genotype,
                FLD_VCF_SAMPLE_ID, b -> b.vcfSampleId = alternateValueSource.VcfSampleId
        );
        nameToAlternateIndexInitializer = Map.of(
                "Chromosome", b ->
                {
                    b.chromosome = alternateValueSource.Chromosome;
                },
                "Position", b ->
                {
                    b.position = alternateValueSource.Position;
                },
                "Ref", b -> b.ref = alternateValueSource.Ref
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
