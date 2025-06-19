package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.DiffFunctions.FILTER_DIFF;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
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
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Test;

public class SomaticVariantDataTest
{
    @Test
    public void fullyMatchesSelf()
    {
        var victim = TestSomaticVariantDataBuilder.create();
        var diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(victim));
        assertNull(victim.findMismatch(victim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(victim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, victim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(victim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(victim, MatchLevel.REPORTABLE, diffThresholds, true));
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        var victim = TestSomaticVariantDataBuilder.create(b ->
        {
            b.comparisonChromosome = "8";
            b.comparisonPosition = 10000;
        });
        var liftoverVictim = TestSomaticVariantDataBuilder.create(b ->
        {
            b.chromosome = "8";
            b.position = 10000;
        });
        var diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, true));
    }

    @Test
    public void onlyMatchesIndex()
    {
        var refVictim = TestSomaticVariantDataBuilder.create();
        var newVictim = TestSomaticVariantDataBuilder.createWithAlternateDefaults(b ->
        {
            b.chromosome = refVictim.Chromosome;
            b.position = refVictim.Position;
            b.ref = refVictim.Ref;
            b.alt = refVictim.Alt;
            b.type = refVictim.Type;
            b.comparisonChromosome = refVictim.mComparisonChromosome;
            b.comparisonPosition = refVictim.mComparisonPosition;
        });

        var diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        var mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(18, mismatch.DiffValues().size());
    }

    @Test
    public void nonPurpleMatchHandledCorrectly()
    {
        var refVictim = TestSomaticVariantDataBuilder.create(b -> b.hasPurpleAnnotation = false);
        var newVictim = TestSomaticVariantDataBuilder.createWithAlternateDefaults(b ->
        {
            b.chromosome = refVictim.Chromosome;
            b.position = refVictim.Position;
            b.ref = refVictim.Ref;
            b.alt = refVictim.Alt;
            b.type = refVictim.Type;
            b.comparisonChromosome = refVictim.mComparisonChromosome;
            b.comparisonPosition = refVictim.mComparisonPosition;
            b.hasPurpleAnnotation = false;
        });

        var diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        var mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(12, mismatch.DiffValues().size());
    }

    @Test
    public void unfilteredMatchHandledCorrectly()
    {
        var passVictim = TestSomaticVariantDataBuilder.create();
        var filteredVictim = TestSomaticVariantDataBuilder.createWithAlternateDefaults(b ->
        {
            b.chromosome = passVictim.Chromosome;
            b.position = passVictim.Position;
            b.ref = passVictim.Ref;
            b.alt = passVictim.Alt;
            b.type = passVictim.Type;
            b.comparisonChromosome = passVictim.mComparisonChromosome;
            b.comparisonPosition = passVictim.mComparisonPosition;
            b.isFromUnfilteredVcf = true;
            b.hasPurpleAnnotation = false;
        });

        var diffThresholds = createDefaultThresholds();

        assertTrue(passVictim.matches(filteredVictim));
        var mismatch = passVictim.findMismatch(filteredVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.REF_ONLY, mismatch.MismatchType());
        assertEquals(passVictim, mismatch.RefItem());
        assertEquals(filteredVictim, mismatch.NewItem());
        assertEquals(7, mismatch.DiffValues().size());

        assertTrue(filteredVictim.matches(passVictim));
        var oppositeMismatch = filteredVictim.findMismatch(passVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.NEW_ONLY, oppositeMismatch.MismatchType());
        assertEquals(filteredVictim, oppositeMismatch.RefItem());
        assertEquals(passVictim, oppositeMismatch.NewItem());
        assertEquals(7, oppositeMismatch.DiffValues().size());
    }

    @Test
    public void reportabilityLossHandledCorrectly()
    {
        var passVictim = TestSomaticVariantDataBuilder.create();
        var nonReportableVictim = TestSomaticVariantDataBuilder.createWithAlternateDefaults(b ->
        {
            b.chromosome = passVictim.Chromosome;
            b.position = passVictim.Position;
            b.ref = passVictim.Ref;
            b.alt = passVictim.Alt;
            b.type = passVictim.Type;
            b.comparisonChromosome = passVictim.mComparisonChromosome;
            b.comparisonPosition = passVictim.mComparisonPosition;
            b.reported = false;
        });

        var diffThresholds = createDefaultThresholds();

        assertTrue(passVictim.matches(nonReportableVictim));
        var mismatch = passVictim.findMismatch(nonReportableVictim, MatchLevel.REPORTABLE, diffThresholds, false);

        assertEquals(MismatchType.REF_ONLY, mismatch.MismatchType());
        assertEquals(passVictim, mismatch.RefItem());
        assertEquals(nonReportableVictim, mismatch.NewItem());
        assertEquals(18, mismatch.DiffValues().size());

        assertTrue(nonReportableVictim.matches(passVictim));
        var oppositeMismatch = nonReportableVictim.findMismatch(passVictim, MatchLevel.REPORTABLE, diffThresholds, false);

        assertEquals(MismatchType.NEW_ONLY, oppositeMismatch.MismatchType());
        assertEquals(nonReportableVictim, oppositeMismatch.RefItem());
        assertEquals(passVictim, oppositeMismatch.NewItem());
        assertEquals(18, oppositeMismatch.DiffValues().size());
    }

    @Test
    public void indexMismatchesAreRecognized()
    {
        var victim = TestSomaticVariantDataBuilder.create();
        var alternateVictim = TestSomaticVariantDataBuilder.createWithAlternateDefaults();

        var chromosomeMismatch = TestSomaticVariantDataBuilder.create(b ->
        {
            b.chromosome = alternateVictim.Chromosome;
            b.comparisonChromosome = alternateVictim.mComparisonChromosome;
        });

        assertFalse(victim.matches(chromosomeMismatch));
        assertFalse(chromosomeMismatch.matches(victim));

        var positionMismatch = TestSomaticVariantDataBuilder.create(b ->
        {
            b.position = alternateVictim.Position;
            b.comparisonPosition = alternateVictim.mComparisonPosition;
        });
        assertFalse(victim.matches(positionMismatch));
        assertFalse(positionMismatch.matches(victim));

        var refMismatch = TestSomaticVariantDataBuilder.create(b -> b.ref = alternateVictim.Ref);

        assertFalse(victim.matches(refMismatch));
        assertFalse(refMismatch.matches(victim));

        var altMismatch = TestSomaticVariantDataBuilder.create(b -> b.alt = alternateVictim.Alt);

        assertFalse(victim.matches(altMismatch));
        assertFalse(altMismatch.matches(victim));

        var variantTypeMismatch = TestSomaticVariantDataBuilder.create(b -> b.type = alternateVictim.Type);

        assertFalse(victim.matches(variantTypeMismatch));
        assertFalse(variantTypeMismatch.matches(victim));
    }

    @Test
    public void singleFieldMismatchIsRecognized()
    {
        var alternateValueSource = TestSomaticVariantDataBuilder.createWithAlternateDefaults();

        testSingleFieldMismatch(FLD_QUAL, b -> b.qual = alternateValueSource.Qual);
        testSingleFieldMismatch(FLD_TIER, b -> b.tier = alternateValueSource.Tier);
        testSingleFieldMismatch(FLD_TUMOR_SUPPORTING_READ_COUNT, b -> b.tumorSupportingReadCount =
                alternateValueSource.TumorSupportingReadCount);
        testSingleFieldMismatch(FLD_TUMOR_TOTAL_READ_COUNT, b -> b.tumorTotalReadCount = alternateValueSource.TumorTotalReadCount);
        testSingleFieldMismatch(FLD_GENE, b -> b.gene = alternateValueSource.Gene);
        testSingleFieldMismatch(FLD_CANON_EFFECT, b -> b.canonicalEffect = alternateValueSource.CanonicalEffect);
        testSingleFieldMismatch(FLD_CODING_EFFECT, b -> b.canonicalCodingEffect = alternateValueSource.CanonicalCodingEffect);
        testSingleFieldMismatch(FLD_HGVS_CODING, b -> b.canonicalHgvsCodingImpact = alternateValueSource.CanonicalHgvsCodingImpact);
        testSingleFieldMismatch(FLD_HGVS_PROTEIN, b -> b.canonicalHgvsProteinImpact = alternateValueSource.CanonicalHgvsProteinImpact);
        testSingleFieldMismatch(FLD_HOTSPOT, b -> b.hotspotStatus = alternateValueSource.HotspotStatus);
        testSingleFieldMismatch(FLD_BIALLELIC, b -> b.biallelic = alternateValueSource.Biallelic);
        testSingleFieldMismatch(FLD_OTHER_REPORTED, b -> b.otherReportedEffects = alternateValueSource.OtherReportedEffects);
        testSingleFieldMismatch(FLD_SUBCLONAL_LIKELIHOOD, b -> b.subclonalLikelihood = alternateValueSource.SubclonalLikelihood);
        testSingleFieldMismatch(FLD_VARIANT_COPY_NUMBER, b -> b.variantCopyNumber = alternateValueSource.VariantCopyNumber);
        testSingleFieldMismatch(FLD_PURITY_ADJUSTED_VAF, b -> b.purityAdjustedVaf = alternateValueSource.PurityAdjustedVaf);
        testSingleFieldMismatch(FLD_LPS, b -> b.hasLPS = alternateValueSource.HasLPS);
        testSingleFieldMismatch(FILTER_DIFF, b -> b.filters = alternateValueSource.Filters);
    }

    public void testSingleFieldMismatch(final String field, final Consumer<TestSomaticVariantDataBuilder> initializer)
    {
        var refVictim = TestSomaticVariantDataBuilder.create();
        var newVictim = TestSomaticVariantDataBuilder.create(initializer);

        assertTrue(refVictim.matches(newVictim));
        assertTrue(newVictim.matches(refVictim));

        var diffThresholds = createDefaultThresholds();
        var detailedMismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, detailedMismatch.MismatchType());
        assertEquals(refVictim, detailedMismatch.RefItem());
        assertEquals(newVictim, detailedMismatch.NewItem());
        assertEquals(1, detailedMismatch.DiffValues().size());
        assertEquals(field, detailedMismatch.DiffValues().get(0).split("\\(")[0]);

        var reportedMismatch = refVictim.findMismatch(newVictim, MatchLevel.REPORTABLE, diffThresholds, false);

        assertEquals(MismatchType.VALUE, reportedMismatch.MismatchType());
        assertEquals(refVictim, reportedMismatch.RefItem());
        assertEquals(newVictim, reportedMismatch.NewItem());
        assertEquals(1, reportedMismatch.DiffValues().size());
        assertEquals(field, reportedMismatch.DiffValues().get(0).split("\\(")[0]);
    }

    private static DiffThresholds createDefaultThresholds()
    {
        var config = new ComparConfig();
        var comparer = new SomaticVariantComparer(config);

        var diffThresholds = new DiffThresholds();
        comparer.registerThresholds(diffThresholds);
        return diffThresholds;
    }
}
