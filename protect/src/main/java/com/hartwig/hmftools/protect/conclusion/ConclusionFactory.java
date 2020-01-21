package com.hartwig.hmftools.protect.conclusion;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionFactory {

    private ConclusionFactory() {

    }

    public static String createConclusion(@NotNull String patientPrimaryTumorLocation, double tumorMTL, double tumorMTB, double tumorMSI,
            double chordScore, @NotNull List<ReportableGeneFusion> geneFusions, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull List<? extends Variant> passSomaticVariants, @NotNull List<TemplateConclusion> templateConclusionList, double purity) {
        String conclusion = Strings.EMPTY;
        String enter = " <enter> ";

        String textTumorLocation = createTumorLocationSentense(patientPrimaryTumorLocation);
        conclusion = textTumorLocation + enter;

        for (TemplateConclusion templateConclusion : templateConclusionList) {
            if (TumorMutationalStatus.fromLoad(tumorMTL).equals(TumorMutationalStatus.HIGH)) {
                if (templateConclusion.abberrationGeneSummary().equals("High TMB")) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion += highMTL + enter;
                }
            } else if (MicrosatelliteStatus.fromIndelsPerMb(tumorMSI).equals(MicrosatelliteStatus.MSI)) {
                if (templateConclusion.abberrationGeneSummary().equals("MSI")) {
                    String highMSI = sentenseHighMSI(tumorMSI, templateConclusion);
                    conclusion += highMSI + enter;
                }
            } else if (chordScore >= 0.5) { //TODO create enum
                if (templateConclusion.abberrationGeneSummary().equals("HR-deficient")) {
                    String hrDeficient = sentenseHrDeficient(chordScore, templateConclusion);
                    conclusion += hrDeficient + enter;
                }
            } else if (purity < 0.20) {
                if (templateConclusion.abberrationGeneSummary().equals("LOW PURITY")) {
                    String lowPurity = sentenceLowPurity(purity, templateConclusion);
                    conclusion += lowPurity + enter;
                }
            }
        }

        conclusion += "\n";

        return conclusion;
    }

    private static String createTumorLocationSentense(@NotNull String patientPrimaryTumorLocation) {
        String tumorLocation = Strings.EMPTY;
        return tumorLocation + "sample showing";
    }

    private static String sentenseHighMTL(double tumorMTL, double tumorMTB, @NotNull TemplateConclusion templateConclusion) {
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("mutational load (ML) of XXX", "mutational load (ML) of " + Double.toString(tumorMTL));
        sentence = sentence.replace("tumor mutation burden (TMB) of XXX mt/Mb",
                "tumor mutation burden (TMB) of " + String.format("%.1f", tumorMTB) + " mt/Mb");
        return sentence;

    }

    private static String sentenseHighMSI(double tumorMSI, @NotNull TemplateConclusion templateConclusion) {
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("MSI signature score XXX ", "MSI signature score " + String.format("%.2f", tumorMSI));
        return sentence;
    }

    private static String sentenseHrDeficient(double chordScore, @NotNull TemplateConclusion templateConclusion) {
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("a high CHORD score (XXX) ", "a high CHORD score (" + String.format("%.2f", chordScore) + ") ");
        return sentence;
    }

    private static String sentenceLowPurity(double purity, @NotNull TemplateConclusion templateConclusion) {
        String purityPercentage = new DecimalFormat("#'%'").format(purity * 100);
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("Due to the low tumor purity (XXX) ", "Due to the low tumor purity (" + purityPercentage + ")");
        return sentence;
    }
}
