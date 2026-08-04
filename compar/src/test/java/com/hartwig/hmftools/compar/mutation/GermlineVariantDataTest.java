package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_FILTER;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_QUAL;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_CANON_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_CODING_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_GENE;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_HGVS_CODING;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_HGVS_PROTEIN;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_HOTSPOT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_OTHER_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_PURITY_ADJUSTED_VAF;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_TIER;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_TUMOR_SUPPORTING_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_TUMOR_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_VARIANT_COPY_NUMBER;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.FieldCheckCache;
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
        comparer = new GermlineVariantComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestGermlineVariantDataBuilder.BUILDER;
        GermlineVariantData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FLD_HOTSPOT, b -> b.hotspotStatus = alternateValueSource.HotspotStatus);
        fieldToAlternateValueInitializer.put(FLD_TIER, b -> b.tier = alternateValueSource.Tier);
        fieldToAlternateValueInitializer.put(FLD_GENE, b -> b.gene = alternateValueSource.Gene);
        fieldToAlternateValueInitializer.put(FLD_CANON_EFFECT, b -> b.canonicalEffect = alternateValueSource.CanonicalEffect);
        fieldToAlternateValueInitializer.put(FLD_CODING_EFFECT, b -> b.canonicalCodingEffect =
                alternateValueSource.CanonicalCodingEffect);
        fieldToAlternateValueInitializer.put(FLD_HGVS_CODING, b -> b.canonicalHgvsCodingImpact =
                alternateValueSource.CanonicalHgvsCodingImpact);
        fieldToAlternateValueInitializer.put(FLD_HGVS_PROTEIN, b -> b.canonicalHgvsProteinImpact =
                alternateValueSource.CanonicalHgvsProteinImpact);
        fieldToAlternateValueInitializer.put(FLD_OTHER_REPORTED, b -> b.otherReportedEffects =
                alternateValueSource.OtherReportedEffects);
        fieldToAlternateValueInitializer.put(FLD_QUAL, b -> b.qual = alternateValueSource.Qual);
        fieldToAlternateValueInitializer.put(FLD_VARIANT_COPY_NUMBER, b -> b.variantCopyNumber =
                alternateValueSource.VariantCopyNumber);
        fieldToAlternateValueInitializer.put(FLD_PURITY_ADJUSTED_VAF, b -> b.purityAdjustedVaf =
                alternateValueSource.PurityAdjustedVaf);
        fieldToAlternateValueInitializer.put(FLD_TUMOR_SUPPORTING_READ_COUNT, b -> b.tumorSupportingReadCount =
                alternateValueSource.TumorSupportingReadCount);
        fieldToAlternateValueInitializer.put(FLD_TUMOR_TOTAL_READ_COUNT, b -> b.tumorTotalReadCount =
                alternateValueSource.TumorTotalReadCount);
        fieldToAlternateValueInitializer.put(FLD_FILTER, b -> b.filters = alternateValueSource.Filters);

        nameToAlternateIndexInitializer = Map.of(
                "Chromosome", b ->
                {
                    b.chromosome = alternateValueSource.Chromosome;
                    b.comparisonChromosome = alternateValueSource.mComparisonChromosome;
                },
                "Position", b ->
                {
                    b.position = alternateValueSource.Position;
                    b.comparisonPosition = alternateValueSource.mComparisonPosition;
                },
                "Ref", b -> b.ref = alternateValueSource.Ref,
                "Alt", b -> b.alt = alternateValueSource.Alt,
                "Type", b -> b.type = alternateValueSource.Type
        );
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
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
        FieldCheckCache detailedFieldConfig = createDefaultThresholds(MatchLevel.DETAILED);
        FieldCheckCache reportableFieldConfig = createDefaultThresholds(MatchLevel.REPORTABLE);

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(comparer, liftoverVictim, MatchLevel.DETAILED, false));
        assertNull(victim.findMismatch(comparer, liftoverVictim, MatchLevel.REPORTABLE, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(comparer, liftoverVictim, MatchLevel.DETAILED, true));
        assertEquals(expectedMatch, victim.findMismatch(comparer, liftoverVictim, MatchLevel.REPORTABLE, true));
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
