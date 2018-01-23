package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.report.Commons;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public class ExplanationPage {

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
                        "Findings in the tumor that also exist in the germline are not included in this report.",
                        "Analysis for samples with an implied tumor purity below 20% suffer from lower sensitivity. "
                                + "This means we have a lower likelihood of finding real variants and "
                                + "we likely underestimate the true mutational load."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> snvIndelExplanationSection() {
        return toList("Details on reported genomic variant fields",
                Lists.newArrayList(
                        "The 'position' refers to the chromosome and start base of the variant with " + "respect to this reference genome.",
                        "The 'variant' displays what was expected as reference base and what " + "was found instead ('ref' > 'alt').",
                        "The 'depth (VAF)' displays the number of observations of the specific variant versus "
                                + "the total number of reads in this location in the format 'alt / total (%)'.",
                        "The 'predicted effect' provides additional information on the variant, including "
                                + "the change in coding sequence ('c.'), the change in protein ('p.') and "
                                + "the predicted impact on the final protein on the second line of this field.",
                        "The 'cosmic' fields display a link to the COSMIC database which contains "
                                + "additional information on the variant. If the variant could not be found in the "
                                + "COSMIC database, this field will be left blank. The CosmicGenes v76 database is used "
                                + "to look-up these IDs.",
                        "The implied tumor purity is the percentage of tumor DNA in the biopsy based on analysis of "
                                + "whole genome data.",
                        "The 'Ploidy (TAF)' field displays the tumor ploidy for the observed variant. The ploidy "
                                + "has been adjusted for the implied tumor purity (see above) and is shown as a "
                                + "proportion of A’s and B’s (e.g. AAABB for 3 copies A, and 2 copies B). "
                                + "The copy number is the sum of A’s and B’s. The TAF (Tumor adjusted Alternative "
                                + "Frequency) value refers to the alternative allele frequency after correction " + "for tumor purity.",
                        "The tumor mutational load is the total number of somatic missense variants found across"
                                + " the whole genome of the tumor biopsy."));
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
