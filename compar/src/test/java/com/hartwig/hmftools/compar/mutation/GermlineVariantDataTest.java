package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.FILTER_DIFF;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_BIALLELIC;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CANON_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CODING_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_GENE;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_CODING;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_PROTEIN;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HOTSPOT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_OTHER_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_PURITY_ADJUSTED_VAF;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TIER;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_SUPPORTING_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_VARIANT_COPY_NUMBER;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class GermlineVariantDataTest extends ComparableItemTest<GermlineVariantData, GermlineVariantComparer, TestGermlineVariantDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GermlineVariantComparer(new ComparConfig());
        builder = TestGermlineVariantDataBuilder.BUILDER;
        GermlineVariantData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FLD_HOTSPOT, b -> b.hotspotStatus = alternateValueSource.Variant.hotspot());
        fieldToAlternateValueInitializer.put(FLD_TIER, b -> b.tier = alternateValueSource.Variant.tier());
        fieldToAlternateValueInitializer.put(FLD_BIALLELIC, b -> b.biallelic = alternateValueSource.Variant.biallelic());
        fieldToAlternateValueInitializer.put(FLD_GENE, b -> b.gene = alternateValueSource.Variant.gene());
        fieldToAlternateValueInitializer.put(FLD_CANON_EFFECT, b -> b.canonicalEffect = alternateValueSource.Variant.canonicalEffect());
        fieldToAlternateValueInitializer.put(FLD_CODING_EFFECT, b -> b.canonicalCodingEffect =
                alternateValueSource.Variant.canonicalCodingEffect());
        fieldToAlternateValueInitializer.put(FLD_HGVS_CODING, b -> b.canonicalHgvsCodingImpact =
                alternateValueSource.Variant.canonicalHgvsCodingImpact());
        fieldToAlternateValueInitializer.put(FLD_HGVS_PROTEIN, b -> b.canonicalHgvsProteinImpact =
                alternateValueSource.Variant.canonicalHgvsProteinImpact());
        fieldToAlternateValueInitializer.put(FLD_OTHER_REPORTED, b -> b.otherReportedEffects =
                alternateValueSource.Variant.otherReportedEffects());
        fieldToAlternateValueInitializer.put(FLD_QUAL, b -> b.qual = alternateValueSource.Variant.qual());
        fieldToAlternateValueInitializer.put(FLD_VARIANT_COPY_NUMBER, b -> b.variantCopyNumber =
                alternateValueSource.Variant.variantCopyNumber());
        fieldToAlternateValueInitializer.put(FLD_PURITY_ADJUSTED_VAF, b -> b.purityAdjustedVaf =
                alternateValueSource.Variant.adjustedVAF());
        fieldToAlternateValueInitializer.put(FLD_TUMOR_SUPPORTING_READ_COUNT, b -> b.tumorSupportingReadCount =
                alternateValueSource.Variant.allelicDepth().AlleleReadCount);
        fieldToAlternateValueInitializer.put(FLD_TUMOR_TOTAL_READ_COUNT, b -> b.tumorTotalReadCount =
                alternateValueSource.Variant.allelicDepth().TotalReadCount);
        fieldToAlternateValueInitializer.put(FILTER_DIFF, b -> b.filters = alternateValueSource.Filters);

        nameToAlternateIndexInitializer = Map.of(
                "Chromosome", b ->
                {
                    b.chromosome = alternateValueSource.Variant.chromosome();
                    b.comparisonChromosome = alternateValueSource.mComparisonPosition.Chromosome;
                },
                "Position", b ->
                {
                    b.position = alternateValueSource.Variant.position();
                    b.comparisonPosition = alternateValueSource.mComparisonPosition.Position;
                },
                "Ref", b -> b.ref = alternateValueSource.Variant.ref(),
                "Alt", b -> b.alt = alternateValueSource.Variant.alt(),
                "Type", b -> b.type = alternateValueSource.Variant.type()
        );
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        GermlineVariantData victim = TestGermlineVariantDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosome = "8";
            b.comparisonPosition = 10000;
        });
        GermlineVariantData liftoverVictim = TestGermlineVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "8";
            b.position = 10000;
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
        GermlineVariantData victim = TestGermlineVariantDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosome = "8";
            b.comparisonPosition = 10000;
        });
        assertFalse(victim.key().isEmpty());
    }
}
