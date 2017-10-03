package com.hartwig.hmftools.strelka;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.slicing.Slicer;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

final class StrelkaPostProcess {
    private static final Logger LOGGER = LogManager.getLogger(StrelkaPostProcess.class);

    private static final String SNP_QUAL_FIELD = "QSS_NT";
    private static final String INDEL_QUAL_FIELD = "QSI_NT";
    private static final String SNP_TIER_INDEX_FIELD = "TQSS_NT";
    private static final String INDEL_TIER_INDEX_FIELD = "TQSI_NT";
    private static final String TIR_FIELD = "TIR";
    private static final String TAR_FIELD = "TAR";
    private static final String SEPARATOR = ",";
    private static final String TUMOR_GENOTYPE = "TUMOR";

    private static final double THRESHOLD = 1.3;
    private static final int LC_QUALITY_SCORE_THRESHOLD = 20;
    private static final double LC_ALLELE_FREQUENCY_THRESHOLD = 0.1;

    private StrelkaPostProcess() {
    }

    @VisibleForTesting
    static boolean checkVariant(@NotNull final VariantContext variant, @NotNull final Slicer highConfidenceSlicer) {
        if (variant.getAlleles().size() > 2) {
            LOGGER.warn("More than 1 alt for record {}: {}", variant.getContig(), variant.getStart());
            return true;
        } else if (variant.getAlleles().size() < 2) {
            LOGGER.warn("Alt is . for record {}: {}", variant.getContig(), variant.getStart());
            return false;
        }
        try {
            return qualityScore(variant) > LC_QUALITY_SCORE_THRESHOLD && allelicFrequency(variant) > LC_ALLELE_FREQUENCY_THRESHOLD || (
                    highConfidenceSlicer.includes(variantGenomePosition(variant))
                            && allelicFrequency(variant) * qualityScore(variant) > THRESHOLD);
        } catch (final HartwigException e) {
            LOGGER.error("encountered error while processing variant {}: {}:\t{}", variant.getContig(), variant.getStart(), e.getMessage());
            return false;
        }
    }

    @VisibleForTesting
    static int qualityScore(@NotNull final VariantContext variant) throws HartwigException {
        if (variant.isSNP()) {
            return getIntField(variant, SNP_QUAL_FIELD);
        } else if (variant.isIndel()) {
            return getIntField(variant, INDEL_QUAL_FIELD);
        } else {
            throw new HartwigException("record is not indel or snp: " + variant);
        }
    }

    private static int getIntField(@NotNull final VariantContext variant, @NotNull final String fieldKey) {
        try {
            return Integer.parseInt(variant.getAttributeAsString(fieldKey, ""));
        } catch (final NumberFormatException e) {
            LOGGER.warn("Could not parse integer attribute {}.", fieldKey);
            throw e;
        }
    }

    @VisibleForTesting
    static double allelicFrequency(@NotNull final VariantContext variant) throws HartwigException {
        final Genotype tumorGenotype = variant.getGenotype(TUMOR_GENOTYPE);
        if (variant.isSNP()) {
            final int tierIndex = getIntField(variant, SNP_TIER_INDEX_FIELD) - 1;
            final String altField = tumorGenotype.getExtendedAttribute(alleleKey(variant.getAlternateAllele(0))).toString();
            final String refField = tumorGenotype.getExtendedAttribute(alleleKey(variant.getReference())).toString();
            return readAf(altField, refField, tierIndex);
        } else if (variant.isIndel()) {
            final int tierIndex = getIntField(variant, INDEL_TIER_INDEX_FIELD) - 1;
            final String altField = tumorGenotype.getExtendedAttribute(TIR_FIELD).toString();
            final String refField = tumorGenotype.getExtendedAttribute(TAR_FIELD).toString();
            return readAf(altField, refField, tierIndex);
        } else {
            throw new HartwigException("record is not indel or snp: " + variant);
        }
    }

    private static double readAf(@NotNull final String altField, @NotNull final String refField, final int tierIndex)
            throws HartwigException {
        final String[] altFieldParts = altField.split(SEPARATOR);
        final String[] refFieldParts = refField.split(SEPARATOR);
        final double altReads = Double.parseDouble(altFieldParts[tierIndex]);
        final double refReads = Double.parseDouble(refFieldParts[tierIndex]);
        final double total = altReads + refReads;
        if (total == 0) {
            return 0;
        }
        return altReads / total;
    }

    static VariantContext simplifyVariant(@NotNull final VariantContext variant, @NotNull final String sampleName) throws HartwigException {
        //MIVO: force GT to 0/1 even for variants with multiple alts
        final List<Allele> outputVariantAlleles = variant.getAlleles().subList(0, 2);
        final Genotype genotype = new GenotypeBuilder(sampleName, outputVariantAlleles).DP(getDP(variant)).AD(getAD(variant)).make();
        return new VariantContextBuilder(variant).genotypes(genotype).make();
    }

    @VisibleForTesting
    static int getDP(@NotNull final VariantContext variant) {
        return variant.getGenotype(TUMOR_GENOTYPE).getDP();
    }

    @VisibleForTesting
    static int[] getAD(@NotNull final VariantContext variant) throws HartwigException {
        final Genotype tumorGenotype = variant.getGenotype(TUMOR_GENOTYPE);
        if (variant.isSNP()) {
            return readAD(variant.getAlleles(), tumorGenotype);
        } else if (variant.isIndel()) {
            final String refADString = tumorGenotype.getExtendedAttribute(TAR_FIELD).toString().split(SEPARATOR)[0];
            final int refAD = Integer.parseInt(refADString);
            final String altADString = tumorGenotype.getExtendedAttribute(TIR_FIELD).toString().split(SEPARATOR)[0];
            final int altAD = Integer.parseInt(altADString);
            return new int[] { refAD, altAD };
        } else {
            throw new HartwigException("record is not indel or snp: " + variant);
        }
    }

    @NotNull
    private static int[] readAD(@NotNull final List<Allele> alleles, @NotNull final Genotype tumorGenotype) {
        final List<Integer> alleleAds = alleles.stream().map(allele -> {
            final String[] alleleADs = tumorGenotype.getExtendedAttribute(alleleKey(allele), "0").toString().split(SEPARATOR);
            if (alleleADs.length > 0) {
                try {
                    return Integer.parseInt(alleleADs[0]);
                } catch (final NumberFormatException e) {
                    return 0;
                }
            }
            return 0;
        }).collect(Collectors.toList());
        final Integer[] ads = new Integer[alleleAds.size()];
        return ArrayUtils.toPrimitive(alleleAds.toArray(ads));
    }

    private static String alleleKey(@NotNull final Allele allele) {
        return allele.getBaseString() + "U";
    }

    @VisibleForTesting
    static GenomePosition variantGenomePosition(@NotNull final VariantContext variant) {
        return new GenomePosition() {
            @Override
            @NotNull
            public String chromosome() {
                return variant.getContig();
            }

            @Override
            public long position() {
                return variant.getStart();
            }
        };
    }
}
