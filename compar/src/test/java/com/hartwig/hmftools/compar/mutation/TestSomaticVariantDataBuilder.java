package com.hartwig.hmftools.compar.mutation;

import java.util.Collections;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestSomaticVariantDataBuilder
{
    public String chromosome = "7";
    public int position = 140453136;
    public String ref = "A";
    public String alt = "T";
    public VariantType type = VariantType.SNP;
    public String gene = "BRAF";
    public boolean reported = true;
    public Hotspot hotspotStatus = Hotspot.HOTSPOT;
    public VariantTier tier = VariantTier.HOTSPOT;
    public boolean biallelic = false;
    public String canonicalEffect = "missense_variant";
    public String canonicalCodingEffect = "MISSENSE";
    public String canonicalHgvsCodingImpact = "c.1799T>A";
    public String canonicalHgvsProteinImpact = "p.Val600Glu";
    public String otherReportedEffects = null;
    public boolean hasLPS = false;
    public int qual = 275;
    public double subclonalLikelihood = 0.;
    public Set<String> filters = Collections.emptySet();
    public double variantCopyNumber = 1.1;
    public double purityAdjustedVaf = 0.45;
    public int tumorTotalReadCount = 116;
    public int tumorSupportingReadCount = 21;
    public boolean isFromUnfilteredVcf = false;
    public boolean hasPurpleAnnotation = true;
    public String comparisonChromosome = "7";
    public int comparisonPosition = 140453136;

    private static final Consumer<TestSomaticVariantDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.chromosome = "8";
        b.position = 10000;
        b.ref = "C";
        b.alt = "G";
        b.type = VariantType.INDEL;
        b.gene = "BRCA1";
        b.reported = false;
        b.hotspotStatus = Hotspot.NEAR_HOTSPOT;
        b.tier = VariantTier.PANEL;
        b.biallelic = true;
        b.canonicalEffect = "synonymous_variant";
        b.canonicalCodingEffect = "SYNONYMOUS";
        b.canonicalHgvsCodingImpact = "c.1800T>A";
        b.canonicalHgvsProteinImpact = "p.Val601Glu";
        b.otherReportedEffects = "OTHER_EFFECT";
        b.hasLPS = true;
        b.qual = 512;
        b.subclonalLikelihood = 1.;
        b.filters = Set.of("PON");
        b.variantCopyNumber = 3.6;
        b.purityAdjustedVaf = 1.1;
        b.tumorTotalReadCount = 312;
        b.tumorSupportingReadCount = 50;
        b.comparisonChromosome = "8";
        b.comparisonPosition = 10000;
    };

    public static final TestComparableItemBuilder<TestSomaticVariantDataBuilder, SomaticVariantData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestSomaticVariantDataBuilder::new,
                    TestSomaticVariantDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private SomaticVariantData build()
    {
        return new SomaticVariantData(chromosome, position, ref, alt, type, gene, reported, hotspotStatus, tier, biallelic, canonicalEffect,
                canonicalCodingEffect, canonicalHgvsCodingImpact, canonicalHgvsProteinImpact, otherReportedEffects, hasLPS, qual,
                subclonalLikelihood, filters, variantCopyNumber, purityAdjustedVaf, tumorSupportingReadCount, tumorTotalReadCount,
                isFromUnfilteredVcf, hasPurpleAnnotation, comparisonChromosome, comparisonPosition);
    }
}
