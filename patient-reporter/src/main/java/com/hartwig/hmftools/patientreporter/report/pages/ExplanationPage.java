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
                cmp.text(Commons.TITLE + " - Report Explanation").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                generalExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                snvIndelExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                copyNumberExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                disruptionExplanationSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                fusionExplanation());
    }

    @NotNull
    private static ComponentBuilder<?, ?> generalExplanationSection() {
        return toList("Details on the report in general",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "Somatic (tumor derived) variation is reported. "
                                + "Potential variants in the tumor that also exist in germline findings are not included in this report.",
                        "Variant detection in samples with lower tumor content is less sensitive. "
                                + "In case of a low tumor purity (below 20%) likelihood of failing to detect potential variants increases.",
                        "The implied tumor purity is the percentage of tumor DNA in the biopsy based on analysis of "
                                + "whole genome data."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> snvIndelExplanationSection() {
        return toList("Details on reported point mutations",
                Lists.newArrayList(
                        "The 'Position' refers to the chromosome and start base of the variant with " + "respect to this reference genome.",
                        "The 'Variant' displays what was expected as reference base and what " + "was found instead ('ref' > 'alt').",
                        "The 'Depth (VAF)' displays the number of observations of the specific variant versus "
                                + "the total number of reads in this location in the format 'alt / total (%)'.",
                        "The 'Predicted Effect' provides additional information on the variant, including "
                                + "the change in coding sequence ('c.'), the change in protein ('p.') and "
                                + "the predicted impact on the final protein on the second line of this field.",
                        "The 'Cosmic' fields display a link to the COSMIC database which contains "
                                + "additional information on the variant. If the variant could not be found in the "
                                + "COSMIC database, this field will be left blank. The CosmicGenes v76 database is used "
                                + "to look-up these IDs.",
                        "The 'Ploidy (TAF)' field displays the tumor ploidy for the observed variant. The ploidy "
                                + "has been adjusted for the implied tumor purity (see above) and is shown as a "
                                + "proportion of A’s and B’s (e.g. AAABB for 3 copies A, and 2 copies B). "
                                + "The copy number is the sum of A’s and B’s. The TAF (Tumor adjusted Alternative "
                                + "Frequency) value refers to the alternative allele frequency after correction " + "for tumor purity."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberExplanationSection() {
        return toList("Details on reported gene copy numbers",
                Lists.newArrayList("The lowest copy number value along the exonic regions of the canonical transcript is determined as "
                                + "a measure for the gene's copy number.",
                        "Copy numbers are corrected for the implied tumor purity and represent the number of copies in the tumor DNA.",
                        "Any gene with no copies is reported as loss.",
                        "Any gene with at least 8 copies is reported as a gain.",
                        "Any gene with more copies than 2.2 times the average tumor ploidy is reported as a gain."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> disruptionExplanationSection() {
        return toList("Details on reported gene disruptions",
                Lists.newArrayList("Genes are only reported as disrupted if their canonical transcript has been disrupted"));
    }

    @NotNull
    private static ComponentBuilder<?, ?> fusionExplanation() {
        return toList("Details on reported gene fusions",
                Lists.newArrayList("Only intronic in-frame fusions or whole exon deletions are reported.",
                        "The canonical, or otherwise longest transcript validly fused is reported.",
                        "Fusions are restricted to those in the Fusion Gene list curated by COSMIC.",
                        "We additionally select fusions where one partner occurs in the 5' or 3' position in COSMIC >3 times.",
                        "See http://cancer.sanger.ac.uk/cosmic/fusion for more information."));
    }
}
