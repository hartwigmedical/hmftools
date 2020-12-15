package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_MINOR_ALLELE_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG;

import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextDecorator implements GenomePosition {

    private final VariantContext context;
    private final VariantType type;
    private final String filter;
    private final String ref;
    private final String alt;
    private final VariantTier tier;
    private final SnpEffSummary snpEffSummary;

    public VariantContextDecorator(final VariantContext context) {
        this.context = context;
        this.filter = displayFilter(context);
        this.type = VariantType.type(context);
        this.ref = context.getReference().getBaseString();
        this.alt = context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));
        this.tier = VariantTier.fromContext(context);
        this.snpEffSummary = SnpEffSummaryFactory.fromSnpEffEnrichment(context);
    }

    @NotNull
    public VariantContext context() {
        return context;
    }

    @NotNull
    public String filter() {
        return filter;
    }

    @NotNull
    @Override
    public String chromosome() {
        return context.getContig();
    }

    @Override
    public long position() {
        return context.getStart();
    }

    @NotNull
    public VariantType type() {
        return type;
    }

    @NotNull
    public String ref() {
        return ref;
    }

    @NotNull
    public String alt() {
        return alt;
    }

    @NotNull
    public SnpEffSummary snpEffSummary() {
        return snpEffSummary;
    }

    @NotNull
    public String gene() {
        return snpEffSummary.gene();
    }

    public double qual() {
        return context.getPhredScaledQual();
    }

    public double adjustedCopyNumber() {
        return context.getAttributeAsDouble(PURPLE_CN_INFO, 0);
    }

    public double adjustedVaf() {
        return context.getAttributeAsDouble(PURPLE_AF_INFO, 0);
    }

    public boolean biallelic() {
        return context.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false);
    }

    public double minorAlleleCopyNumber() {
        return context.getAttributeAsDouble(PURPLE_MINOR_ALLELE_CN_INFO, context.getAttributeAsDouble(PURPLE_MINOR_ALLELE_PLOIDY_INFO, 0));
    }

    public double variantCopyNumber() {
        return context.getAttributeAsDouble(PURPLE_VARIANT_CN_INFO, context.getAttributeAsDouble(PURPLE_VARIANT_PLOIDY_INFO, 0));
    }

    @NotNull
    public AllelicDepth allelicDepth(@NotNull final String sample) {
        final Genotype genotype = context.getGenotype(sample);
        return genotype != null ? AllelicDepth.fromGenotype(genotype) : NO_DEPTH;
    }

    @NotNull
    public GenotypeStatus genotypeStatus(@NotNull final String sample) {
        final Genotype genotype = context.getGenotype(sample);
        return genotype != null ? GenotypeStatus.fromGenotype(genotype) : GenotypeStatus.UNKNOWN;
    }

    @NotNull
    public VariantTier tier() {
        return tier;
    }

    @NotNull
    public PathogenicSummary pathogenicSummary() {
        return PathogenicSummaryFactory.fromContext(context);
    }

    public int repeatCount() {
        return context.getAttributeAsInt(REPEAT_COUNT_FLAG, 0);
    }

    @NotNull
    public String repeatSequence() {
        return context.getAttributeAsString(REPEAT_SEQUENCE_FLAG, Strings.EMPTY);
    }

    @NotNull
    public Hotspot hotspot() {
        return HotspotEnrichment.fromVariant(context);
    }

    @NotNull
    public String trinucleotideContext() {
        return context.getAttributeAsString(TRINUCLEOTIDE_FLAG, Strings.EMPTY);
    }

    public double mappability() {
        return context.getAttributeAsDouble(MAPPABILITY_TAG, 0);
    }

    public boolean reported() {
        return context.getAttributeAsBoolean(REPORTED_FLAG, false);
    }

    @NotNull
    public String microhomology() {
        return context.getAttributeAsString(MICROHOMOLOGY_FLAG, Strings.EMPTY);
    }

    private static String displayFilter(@NotNull final VariantContext context) {
        if (context.isFiltered()) {
            StringJoiner joiner = new StringJoiner(";");
            context.getFilters().forEach(joiner::add);
            return joiner.toString();
        } else {
            return SomaticVariantFactory.PASS_FILTER;
        }
    }

    private static final AllelicDepth NO_DEPTH = new AllelicDepth() {
        @Override
        public int totalReadCount() {
            return 0;
        }

        @Override
        public int alleleReadCount() {
            return 0;
        }
    };
}
