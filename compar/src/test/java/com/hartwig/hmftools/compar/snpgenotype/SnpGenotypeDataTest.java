package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeData.FLD_GENOTYPE;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeData.FLD_VCF_SAMPLE_ID;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
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
        comparer = new SnpGenotypeComparer(new ComparConfig());
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
                    b.comparisonChromosome = alternateValueSource.mComparisonPosition.Chromosome;
                },
                "Position", b ->
                {
                    b.position = alternateValueSource.Position;
                    b.comparisonPosition = alternateValueSource.mComparisonPosition.Position;
                },
                "Ref", b -> b.ref = alternateValueSource.Ref
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        SnpGenotypeData victim = TestSnpGenotypeDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosome = "11";
            b.comparisonPosition = 30000;
        });
        SnpGenotypeData liftoverVictim = TestSnpGenotypeDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "11";
            b.position = 30000;
        });
        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, true));
    }

    @Test
    public void keyNonEmptyForLiftover()
    {
        SnpGenotypeData victim = TestSnpGenotypeDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosome = "11";
            b.comparisonPosition = 30000;
        });
        assertFalse(victim.key().isEmpty());
    }
}
