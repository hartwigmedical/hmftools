package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.report.Commons;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class ExplanationPage {

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(Commons.TITLE_SEQUENCE + " - Report Explanation").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                generalExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                evidenceExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                snvIndelExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                copyNumberExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                fusionExplanation(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                disruptionExplanationSection());
    }

    @NotNull
    private static ComponentBuilder<?, ?> generalExplanationSection() {
        return toList("Details on the report in general",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "Transcripts used for reporting can be found on https://github.com/hartwigmedical and "
                                + "are generally the canonical transcripts as defined by Ensembl.",
                        "Variant detection in samples with lower tumor content is less sensitive. "
                                + "In case of a low tumor purity (below 20%) likelihood of failing to detect potential variants increases.",
                        "The implied tumor purity is the percentage of tumor cells in the biopsy based on analysis of "
                                + "whole genome data."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceExplanationSection() {
        return toList("Details on the reported clinical evidence",
                Lists.newArrayList("CGI, OncoKb and CiViC are used to annotate variants of all types with clinical evidence.",
                        "More information on (CGI) biomarkers can be found on https://www.cancergenomeinterpreter.org/biomarkers",
                        "Clinical trials are matched against the iclusion database (https://iclusion.org)"));
    }

    @NotNull
    private static ComponentBuilder<?, ?> snvIndelExplanationSection() {
        return toList("Details on reported somatic variants",
                Lists.newArrayList("The 'Read Depth' displays the raw number of reads supporting the "
                                + "variant versus the total number of reads on the mutated position.",
                        "The 'Ploidy (VAF)' field displays the tumor ploidy for the observed position. The ploidy "
                                + "has been adjusted for the implied tumor purity (see above) and is shown as a "
                                + "proportion of A’s and B’s (e.g. AAABB for 3 copies of A, and 2 copies of B). "
                                + "The copy number is the sum of A’s and B’s. The VAF value is to the alternative allele "
                                + "frequency after correction for tumor purity.",
                        "The 'Biallelic' field indicates whether the variant is present across all alleles in the tumor.",
                        "The 'Driver' field is based on the driver probability calculated based on the HMF database. "
                                + "A variant with High driver likelihood is likely to be positively selected for."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberExplanationSection() {
        return toList("Details on reported gene copy numbers",
                Lists.newArrayList("The lowest copy number value along the exonic regions of the canonical transcript is determined as "
                                + "a measure for the gene's copy number.",
                        "Copy numbers are corrected for the implied tumor purity and represent the number of copies in the tumor DNA.",
                        "Any gene with less than 0.5 copies along the entire canonical transcript is reported as a full loss. ",
                        "Any gene where only a part along the canonical transcript has less than 0.5 copies is reported as a partial loss. ",
                        "Any gene with more copies than 3 times the average tumor ploidy is reported as a gain."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> fusionExplanation() {
        return toList("Details on reported gene fusions",
                Lists.newArrayList("Only intronic in-frame fusions or whole exon deletions are reported.",
                        "The canonical, or otherwise longest transcript validly fused is reported.",
                        "Fusions are restricted to those in a known fusion list based on CiViC, OncoKB, CGI and COSMIC",
                        "We additionally select fusions where one partner is promiscuous in either 5' or 3' position."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> disruptionExplanationSection() {
        return toList("Details on reported gene disruptions",
                Lists.newArrayList("Genes are reported as being disrupted if their canonical transcript has been disrupted",
                        "The range of the disruption is indicated by the intron/exon/promoter region of the break point occurred "
                                + "and the direction the disruption faces.",
                        "The type of disruption can be INV (inversion), DEL (deletion), DUP (duplication), "
                                + "INS (insertion) or BND (translocation)."));
    }
}
