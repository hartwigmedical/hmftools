package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.ComparTestUtil.assertDifferencesAreForFields;
import static com.hartwig.hmftools.compar.ComparTestUtil.union;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_FILTER;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_QUAL;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_BIALLELIC;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_BIALLELIC_PROB;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
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
import java.util.Set;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class SomaticVariantDataTest extends ComparableItemTest<SomaticVariantData, SomaticVariantComparer, TestSomaticVariantDataBuilder>
{
    private static final Set<String> PAVE_ONLY_FIELDS =
            Set.of(FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT, FLD_HGVS_CODING, FLD_HGVS_PROTEIN);
    private static final Set<String> SAGE_ONLY_FIELDS =
            Set.of(FLD_QUAL, FLD_REPORTED, FLD_TIER, FLD_TUMOR_SUPPORTING_READ_COUNT, FLD_TUMOR_TOTAL_READ_COUNT, FLD_LPS, FLD_FILTER);
    private static final Set<String> FIELDS_UP_TO_PAVE = union(SAGE_ONLY_FIELDS, PAVE_ONLY_FIELDS);

    @Before
    public void setUp()
    {
        comparer = new SomaticVariantComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestSomaticVariantDataBuilder.BUILDER;
        SomaticVariantData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FLD_QUAL, b -> b.qual = alternateValueSource.Qual);
        fieldToAlternateValueInitializer.put(FLD_TIER, b -> b.tier = alternateValueSource.Tier);
        fieldToAlternateValueInitializer.put(FLD_TUMOR_SUPPORTING_READ_COUNT, b -> b.tumorSupportingReadCount =
                alternateValueSource.TumorSupportingReadCount);
        fieldToAlternateValueInitializer.put(FLD_TUMOR_TOTAL_READ_COUNT, b -> b.tumorTotalReadCount =
                alternateValueSource.TumorTotalReadCount);
        fieldToAlternateValueInitializer.put(FLD_GENE, b -> b.gene = alternateValueSource.Gene);
        fieldToAlternateValueInitializer.put(FLD_CANON_EFFECT, b -> b.canonicalEffect = alternateValueSource.CanonicalEffect);
        fieldToAlternateValueInitializer.put(FLD_CODING_EFFECT, b -> b.canonicalCodingEffect = alternateValueSource.CanonicalCodingEffect);
        fieldToAlternateValueInitializer.put(FLD_HGVS_CODING, b -> b.canonicalHgvsCodingImpact =
                alternateValueSource.CanonicalHgvsCodingImpact);
        fieldToAlternateValueInitializer.put(FLD_HGVS_PROTEIN, b -> b.canonicalHgvsProteinImpact =
                alternateValueSource.CanonicalHgvsProteinImpact);
        fieldToAlternateValueInitializer.put(FLD_HOTSPOT, b -> b.hotspotStatus = alternateValueSource.HotspotStatus);
        fieldToAlternateValueInitializer.put(FLD_BIALLELIC, b -> b.biallelic = alternateValueSource.Biallelic);
        fieldToAlternateValueInitializer.put(FLD_BIALLELIC_PROB, b -> b.biallelicProb = alternateValueSource.BiallelicProbability);
        fieldToAlternateValueInitializer.put(FLD_OTHER_REPORTED, b -> b.otherReportedEffects = alternateValueSource.OtherReportedEffects);
        fieldToAlternateValueInitializer.put(FLD_SUBCLONAL_LIKELIHOOD, b -> b.subclonalLikelihood =
                alternateValueSource.SubclonalLikelihood);
        fieldToAlternateValueInitializer.put(FLD_VARIANT_COPY_NUMBER, b -> b.variantCopyNumber = alternateValueSource.VariantCopyNumber);
        fieldToAlternateValueInitializer.put(FLD_PURITY_ADJUSTED_VAF, b -> b.purityAdjustedVaf = alternateValueSource.PurityAdjustedVaf);
        fieldToAlternateValueInitializer.put(FLD_LPS, b -> b.hasLPS = alternateValueSource.HasLPS);
        fieldToAlternateValueInitializer.put(FLD_FILTER, b -> b.filters = alternateValueSource.Filters);

        nameToAlternateIndexInitializer = Map.of(
                "Chromosome", b -> b.chromosome = alternateValueSource.Chromosome,
                "Position", b -> b.position = alternateValueSource.Position,
                "Ref", b -> b.ref = alternateValueSource.Ref,
                "Alt", b -> b.alt = alternateValueSource.Alt,
                "Type", b -> b.type = alternateValueSource.Type
        );
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
        nameToNonPassInitializer = Map.of(
                "UnfilteredVcf", b -> b.isFromUnfilteredVcf = true,
                "NonGene", b -> {
                    b.gene = "";
                    b.reported = false;
                }
        );
    }

    @Ignore
    @Test
    public void unfilteredMatchHandledCorrectly()
    {
        SomaticVariantData passVictim = TestSomaticVariantDataBuilder.BUILDER.create();
        SomaticVariantData filteredVictim = TestSomaticVariantDataBuilder.BUILDER.createWithAlternateDefaults(b ->
        {
            b.chromosome = passVictim.Chromosome;
            b.position = passVictim.Position;
            b.ref = passVictim.Ref;
            b.alt = passVictim.Alt;
            b.type = passVictim.Type;
            b.isFromUnfilteredVcf = true;
            b.filters = Set.of("TumorQual");
        });

        assertTrue(passVictim.matches(filteredVictim));

        MatchLevel matchLevel = MatchLevel.DETAILED;
        FieldCheckCache fieldConfig = createDefaultThresholds(matchLevel);
        Mismatch mismatch = passVictim.findMismatch(comparer, filteredVictim, matchLevel, false);

        assertEquals(MismatchType.OLD_ONLY, mismatch.Type);
        assertEquals(passVictim, mismatch.OldItem);
        assertEquals(filteredVictim, mismatch.NewItem);
        assertDifferencesAreForFields(union(SAGE_ONLY_FIELDS, Set.of(FLD_FILTER)), mismatch.DiffValues);

        assertTrue(filteredVictim.matches(passVictim));
        Mismatch oppositeMismatch = filteredVictim.findMismatch(comparer, passVictim, matchLevel, false);

        assertEquals(MismatchType.NEW_ONLY, oppositeMismatch.Type);
        assertEquals(filteredVictim, oppositeMismatch.OldItem);
        assertEquals(passVictim, oppositeMismatch.NewItem);
        assertDifferencesAreForFields(union(SAGE_ONLY_FIELDS, Set.of(FLD_FILTER)), mismatch.DiffValues);
    }
}
