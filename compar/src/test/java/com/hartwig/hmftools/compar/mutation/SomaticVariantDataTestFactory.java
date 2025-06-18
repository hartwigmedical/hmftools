package com.hartwig.hmftools.compar.mutation;

import java.util.Set;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

public record SomaticVariantDataTestFactory(SomaticVariantData variantData)
{
    public static SomaticVariantDataTestFactory createDefault()
    {
        final SomaticVariantData defaultVariantData =
                new SomaticVariantData("7", 140453136, "A", "T", VariantType.SNP, "BRAF", true,
                        Hotspot.HOTSPOT, VariantTier.HOTSPOT, false, "missense_variant",
                        "MISSENSE", "c.1799T>A", "p.Val600Glu",
                        null, false, 275, 0., Set.of("PASS"), 1.1,
                        0.45, new AllelicDepth(116, 21), false,
                        true, "7", 140453136);
        return new SomaticVariantDataTestFactory(defaultVariantData);
    }

    public static SomaticVariantDataTestFactory createAlternateDefault()
    {
        final SomaticVariantData defaultVariantData =
                new SomaticVariantData("8", 10000, "C", "G", VariantType.INDEL, "BRCA1", false,
                        Hotspot.NEAR_HOTSPOT, VariantTier.PANEL, true, "synonymous_variant",
                        "SYNONYMOUS", "c.1800T>A", "p.Val601Glu",
                        "OTHER_EFFECT", true, 512, 1., Set.of("PON"), 3.6,
                        1.1, new AllelicDepth(312, 50), false,
                        true, "8", 10000);
        return new SomaticVariantDataTestFactory(defaultVariantData);
    }

    public SomaticVariantData build()
    {
        return variantData;
    }

    public SomaticVariantDataTestFactory withChromosome(String Chromosome)
    {
        SomaticVariantData newVariant = new SomaticVariantData(Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withPosition(int Position)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withRef(String Ref)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withAlt(String Alt)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withType(VariantType Type)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withReported(boolean Reported)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withCanonicalHgvsProteinImpact(String CanonicalHgvsProteinImpact)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withIsFromUnfilteredVcf(boolean IsFromUnfilteredVcf)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withHasPurpleAnnotation(boolean HasPurpleAnnotation)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, HasPurpleAnnotation, variantData.mComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withComparisonChromosome(String ComparisonChromosome)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, ComparisonChromosome,
                variantData.mComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }

    public SomaticVariantDataTestFactory withComparisonPosition(int ComparisonPosition)
    {
        SomaticVariantData newVariant = new SomaticVariantData(variantData.Chromosome, variantData.Position, variantData.Ref,
                variantData.Alt, variantData.Type, variantData.Gene, variantData.Reported, variantData.HotspotStatus, variantData.Tier,
                variantData.Biallelic, variantData.CanonicalEffect, variantData.CanonicalCodingEffect, variantData.CanonicalHgvsCodingImpact,
                variantData.CanonicalHgvsProteinImpact, variantData.OtherReportedEffects, variantData.HasLPS, variantData.Qual,
                variantData.SubclonalLikelihood, variantData.Filters, variantData.VariantCopyNumber, variantData.PurityAdjustedVaf,
                variantData.TumorDepth, variantData.IsFromUnfilteredVcf, variantData.HasPurpleAnnotation, variantData.mComparisonChromosome,
                ComparisonPosition);
        return new SomaticVariantDataTestFactory(newVariant);
    }
}
