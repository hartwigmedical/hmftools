package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;
import com.hartwig.hmftools.svannotation.GeneAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.TranscriptAnnotation;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class PatientDataSource {

    private static final String COSMIC_IDENTIFIER = "COSM";

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);

    static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    static final FieldBuilder<?> HGVS_CODING_FIELD = field("hgvs_coding", String.class);
    static final FieldBuilder<?> HGVS_PROTEIN_FIELD = field("hgvs_protein", String.class);
    static final FieldBuilder<?> CONSEQUENCE_FIELD = field("consequence", String.class);
    static final FieldBuilder<?> COSMIC_FIELD = field("cosmic", String.class);
    static final FieldBuilder<?> COSMIC_NR_FIELD = field("cosmic_nr", String.class);
    static final FieldBuilder<?> DEPTH_VAF_FIELD = field("depth_vaf", String.class);
    static final FieldBuilder<?> PLOIDY_TAF_FIELD = field("ploidy_taf", String.class);

    static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    static final FieldBuilder<?> BAND_FIELD = field("band", String.class);
    static final FieldBuilder<?> COPY_NUMBER_TYPE_FIELD = field("copynumber_type", String.class);
    static final FieldBuilder<?> COPY_NUMBER_FIELD = field("copynumber", String.class);

    static final FieldBuilder<?> SV_GENE_FIELD = field("gene", String.class);
    static final FieldBuilder<?> SV_POSITION_FIELD = field("breakpoint", String.class);
    static final FieldBuilder<?> SV_TYPE_FIELD = field("type", String.class);
    static final FieldBuilder<?> SV_PARTNER_FIELD = field("partner", String.class);
    static final FieldBuilder<?> SV_HGVS_FIELD = field("hgvs", String.class);
    static final FieldBuilder<?> SV_ORIENTATION_FIELD = field("orientation", String.class);
    static final FieldBuilder<?> SV_GENE_CONTEXT = field("gene context", String.class);
    static final FieldBuilder<?> SV_VAF = field("vaf", String.class);
    static final FieldBuilder<?> SV_TAF = field("taf", String.class);

    private PatientDataSource() {
    }

    @NotNull
    static JRDataSource fromVariants(@NotNull final List<VariantReport> variants, @NotNull final HmfReporterData reporterData) {
        final DRDataSource variantDataSource =
                new DRDataSource(GENE_FIELD.getName(), POSITION_FIELD.getName(), VARIANT_FIELD.getName(), DEPTH_VAF_FIELD.getName(),
                        COSMIC_FIELD.getName(), COSMIC_NR_FIELD.getName(), HGVS_CODING_FIELD.getName(), HGVS_PROTEIN_FIELD.getName(),
                        CONSEQUENCE_FIELD.getName(), PLOIDY_TAF_FIELD.getName());

        for (final VariantReport variant : variants) {
            final String displayGene = reporterData.drupFilter().test(variant) ? variant.gene() + " *" : variant.gene();
            variantDataSource.add(displayGene, variant.chromosomePosition(), variant.variantField(), variant.depthVafField(),
                    variant.cosmicID(), stripCosmicIdentifier(variant.cosmicID()), variant.hgvsCoding(), variant.hgvsProtein(),
                    variant.consequence(), variant.ploidyTafField());
        }

        return variantDataSource;
    }

    @NotNull
    static JRDataSource fromCopyNumbers(@NotNull final List<CopyNumberReport> copyNumbers) {
        final DRDataSource copyNumberDatasource =
                new DRDataSource(CHROMOSOME_FIELD.getName(), BAND_FIELD.getName(), GENE_FIELD.getName(), COPY_NUMBER_TYPE_FIELD.getName(),
                        COPY_NUMBER_FIELD.getName());

        for (final CopyNumberReport copyNumber : copyNumbers) {
            copyNumberDatasource.add(copyNumber.chromosome(), copyNumber.chromosomeBand(), copyNumber.gene(), copyNumber.description(),
                    Integer.toString(copyNumber.copyNumber()));
        }
        return copyNumberDatasource;
    }

    @NotNull
    static JRDataSource fromStructuralVariants(@NotNull List<StructuralVariantAnnotation> variants,
            @NotNull final HmfReporterData reporterData) {

        final DRDataSource svDatasource =
                new DRDataSource(SV_GENE_FIELD.getName(), SV_POSITION_FIELD.getName(), SV_TYPE_FIELD.getName(), SV_PARTNER_FIELD.getName(),
                        SV_HGVS_FIELD.getName(), SV_ORIENTATION_FIELD.getName(), SV_GENE_CONTEXT.getName(), SV_VAF.getName(),
                        SV_TAF.getName());

        final Predicate<GeneAnnotation> inCosmic = g -> reporterData.cosmicModel().data().containsKey(g.getGeneName());
        final Predicate<StructuralVariantAnnotation> intronicDisruption = sv -> {
            for (final GeneAnnotation g : sv.getStartAnnotations().getGenes()) {
                if (sv.getEndAnnotations()
                        .getGenes()
                        .stream()
                        .filter(o -> o.getCanonical().isIntronic() && g.getCanonical().isIntronic()
                                && o.getCanonical().getExonUpstream() == g.getCanonical().getExonUpstream())
                        .count() > 0) {
                    return true;
                }
            }
            return false;
        };

        final List<GeneAnnotation> genes = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : variants) {

            if (intronicDisruption.test(sv)) {
                continue;
            }

            genes.addAll(sv.getStartAnnotations().getGenes());
            genes.addAll(sv.getEndAnnotations().getGenes());
        }

        final ArrayListMultimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        for (final GeneAnnotation g : genes) {
            if (!inCosmic.test(g)) {
                continue;
            }
            geneMap.put(g.getGeneName(), g);
        }

        final Function<TranscriptAnnotation, String> exonDescription = (t) -> {
            if (t.isPromoter()) {
                return "Promoter Region";
            } else if (t.isExonic()) {
                return String.format("In exon %d / %d", t.getExonUpstream(), t.getExonMax());
            } else {
                return String.format("Between exon %d and %d / %d", t.getExonUpstream(), t.getExonDownstream(), t.getExonMax());
            }
        };

        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation g : geneMap.get(geneName)) {
                final StructuralVariant sv = g.getBreakend().getStructuralVariant().getVariant();
                final String partner = sv.startChromosome().equals(sv.endChromosome())
                        ? Long.toString(sv.endPosition() - sv.startPosition())
                        : g.getOtherBreakend().getChromosome();
                final String hgvs = "TODO";
                final String orientation = g.getBreakend().getOrientation() > 0 ? "5\"" : "3\"";
                final String vaf = "TODO";
                final String taf = "TODO";

                svDatasource.add(geneName, g.getBreakend().getPositionString(), sv.type().toString(), partner, hgvs, orientation,
                        exonDescription.apply(g.getCanonical()), vaf, taf);
            }
        }

        return svDatasource;
    }

    @NotNull
    private static String stripCosmicIdentifier(@NotNull final String cosmicID) {
        final int identifierPos = cosmicID.indexOf(COSMIC_IDENTIFIER);
        if (identifierPos >= 0) {
            return cosmicID.substring(identifierPos + COSMIC_IDENTIFIER.length());
        } else {
            return cosmicID;
        }
    }

    @NotNull
    static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, POSITION_FIELD, VARIANT_FIELD, HGVS_CODING_FIELD, HGVS_PROTEIN_FIELD, CONSEQUENCE_FIELD,
                COSMIC_FIELD, COSMIC_NR_FIELD, DEPTH_VAF_FIELD, PLOIDY_TAF_FIELD };
    }

    @NotNull
    static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, BAND_FIELD, GENE_FIELD, COPY_NUMBER_TYPE_FIELD, COPY_NUMBER_FIELD };
    }

    @NotNull
    static FieldBuilder<?>[] structuralVariantFields() {
        return new FieldBuilder<?>[] { SV_GENE_FIELD, SV_POSITION_FIELD, SV_TYPE_FIELD, SV_PARTNER_FIELD, SV_HGVS_FIELD,
                SV_ORIENTATION_FIELD, SV_GENE_CONTEXT, SV_VAF, SV_TAF };
    }
}
