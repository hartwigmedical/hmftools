package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.mutation.GermlineVariantData.FILTER_DELIMITER;

import java.util.Collections;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableGermlineVariantImpl;
import com.hartwig.hmftools.common.variant.ImmutableVariantImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestGermlineVariantDataBuilder
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
    public CodingEffect canonicalCodingEffect = CodingEffect.MISSENSE;
    public String canonicalHgvsCodingImpact = "c.1799T>A";
    public String canonicalHgvsProteinImpact = "p.Val600Glu";
    public String otherReportedEffects = "";
    public double qual = 275;
    public Set<String> filters = Collections.emptySet();
    public double variantCopyNumber = 1.1;
    public double purityAdjustedVaf = 0.45;
    public int tumorTotalReadCount = 116;
    public int tumorSupportingReadCount = 21;
    public String comparisonChromosome = "7";
    public int comparisonPosition = 140453136;

    private static final Consumer<TestGermlineVariantDataBuilder> ALTERNATE_INITIALIZER = b ->
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
        b.canonicalCodingEffect = CodingEffect.SYNONYMOUS;
        b.canonicalHgvsCodingImpact = "c.1800T>A";
        b.canonicalHgvsProteinImpact = "p.Val601Glu";
        b.otherReportedEffects = "OTHER_EFFECT";
        b.qual = 512;
        b.filters = Set.of("minTumorQual");
        b.variantCopyNumber = 3.6;
        b.purityAdjustedVaf = 1.1;
        b.tumorTotalReadCount = 312;
        b.tumorSupportingReadCount = 50;
        b.comparisonChromosome = "8";
        b.comparisonPosition = 10000;
    };

    public static final TestComparableItemBuilder<TestGermlineVariantDataBuilder, GermlineVariantData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGermlineVariantDataBuilder::new,
                    TestGermlineVariantDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GermlineVariantData build()
    {
        GermlineVariant variant = ImmutableGermlineVariantImpl.builder()
                .variant(ImmutableVariantImpl.builder()
                        .chromosome(chromosome)
                        .position(position)
                        .type(type)
                        .gene(gene)
                        .ref(ref)
                        .alt(alt)
                        .allelicDepth(new AllelicDepth(tumorTotalReadCount, tumorSupportingReadCount))
                        .canonicalEffect(canonicalEffect)
                        .canonicalCodingEffect(canonicalCodingEffect)
                        .canonicalHgvsCodingImpact(canonicalHgvsCodingImpact)
                        .canonicalHgvsProteinImpact(canonicalHgvsProteinImpact)
                        .qual(qual)
                        .filter(filters.stream().sorted().collect(Collectors.joining(FILTER_DELIMITER)))
                        .otherReportedEffects(otherReportedEffects)
                        .hotspot(hotspotStatus)
                        .adjustedVAF(purityAdjustedVaf)
                        .variantCopyNumber(variantCopyNumber)
                        .biallelic(biallelic)
                        .reported(reported)
                        .tier(tier)
                        .canonicalTranscript("")
                        .genesAffected(-1)
                        .spliceRegion(false)
                        .worstCodingEffect(CodingEffect.UNDEFINED)
                        .mappability(-1)
                        .adjustedCopyNumber(-1)
                        .minorAlleleCopyNumber(-1)
                        .genotypeStatus(GenotypeStatus.UNKNOWN)
                        .germlineStatus(GermlineStatus.UNKNOWN)
                        .trinucleotideContext("")
                        .microhomology("")
                        .repeatSequence("")
                        .repeatCount(-1)
                        .rnaDepth(null)
                        .clinvarInfo(null)
                        .build())
                .pathogenicity(null)
                .pathogenic(true)
                .build();
        BasePosition comparisonBasePosition = new BasePosition(comparisonChromosome, comparisonPosition);
        return new GermlineVariantData(variant, comparisonBasePosition);
    }
}
