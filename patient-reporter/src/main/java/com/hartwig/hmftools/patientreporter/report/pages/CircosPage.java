package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.PDFWriter.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.ImageBuilder;
import net.sf.dynamicreports.report.builder.component.TextFieldBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalImageAlignment;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class CircosPage {

    @NotNull
    abstract String circosImagePath();

    @NotNull
    private static final String TITLE = "CIRCOS plot";

    @NotNull
    private TextFieldBuilder<String> title() {
        return cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - " + TITLE).setStyle(sectionHeaderStyle());
    }

    @NotNull
    private ImageBuilder circosImage() {
        return cmp.image(circosImagePath()).setHeight(430).setWidth(430).setHorizontalImageAlignment(HorizontalImageAlignment.CENTER);
    }

    public ComponentBuilder<?, ?> reportComponent() {
        // @formatter:off
        return cmp.verticalList(
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    title(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    circosImage(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    description()
                );
        // @formatter:on
    }

    private VerticalListBuilder description() {
        return toList("Details on CIRCOS plot", Lists.newArrayList(
                "The outer first circle shows the chromosomes.  The darker shaded areas represent large "
                        + "gaps in the human reference genome:  i.e. regions of centromeres, heterochromatin & missing short arms.",
                "The second circle shows the allele frequency of all observed somatic variants (incl. exon, intron and intergenic regions). "
                        + "The scale is from 0 to 100%. Each purple dot represents a single somatic variant. "
                        + "Variant allele frequencies have been corrected for tumor purity. "
                        + "Indicated variants are colored according to the type of base change (e.g. C>T/G>A in red) and are in concordance "
                        + "with the coloring used in Alexandrov at al. 2013 Nature paper that describes the use of mutational signatures.",
                "The third circle shows all observed tumor purity adjusted copy number changes, including both focal and chromosomal somatic events. "
                        + "Copy number losses are indicated in red, green shows regions of copy number gain. "
                        + "The scale ranges from 0 (complete loss) to 6 (high level gains). "
                        + "If the absolute copy number is > 6 it is shown as 6 with a green dot on the diagram.",
                "The fourth circle represents the observed 'minor allele copy numbers’ across the chromosome. "
                        + "The range of the chart is from 0 to 3. The expected normal minor allele copy number is 1, "
                        + "and anything below 1 is shown as a loss (orange) and represents a LOH event. "
                        + "Minor allele copy numbers above 1 (blue) indicate amplification events of both A and B alleles at the indicated locations.",
                "The innermost circle displays the observed structural variants within or between the chromosomes. "
                        + "Translocations are indicated in blue, deletions in red, insertions in yellow, tandem duplications in "
                        + "green and inversions in black."));
    }
}
