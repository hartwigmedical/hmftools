package com.hartwig.hmftools.protect.conclusion;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.report.chord.ChordStatus;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionFactory {
    private static final Logger LOGGER = LogManager.getLogger(ConclusionFactory.class);

    private ConclusionFactory() {

    }

    public static StringBuilder createConclusion(@NotNull String patientPrimaryTumorLocation, int tumorMTL, double tumorMTB, double tumorMSI,
            double chordScore, @NotNull List<ReportableGeneFusion> geneFusions, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull List<? extends Variant> passSomaticVariants, @NotNull List<TemplateConclusion> templateConclusionList, double purity,
            @NotNull List<TumorLocationConclusion> tumorLocationConclusion, @NotNull String cancerSubtype) {

        StringBuilder conclusion = new StringBuilder();
        String enter = " <enter> ";

        String textTumorLocation = createTumorLocationSentense(patientPrimaryTumorLocation, tumorLocationConclusion, cancerSubtype);
        conclusion.append(textTumorLocation).append(enter);

        for (TemplateConclusion templateConclusion : templateConclusionList) {
            if (TumorMutationalStatus.fromLoad(tumorMTL).equals(TumorMutationalStatus.HIGH)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.HIGH_MTL)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(highMTL).append(enter);
                }
            } else if (MicrosatelliteStatus.fromIndelsPerMb(tumorMSI).equals(AberrationGenSummary.MSI)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.MSI)) {
                    String highMSI = sentenseHighMSI(tumorMSI, templateConclusion);
                    conclusion.append(highMSI).append(enter);
                }
            } else if (ChordStatus.formChord(chordScore).equals(ChordStatus.HR_DEFICIENT)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary())
                        .equals(AberrationGenSummary.HR_DEFICIENT)) {
                    String hrDeficient = sentenseHrDeficient(chordScore, templateConclusion);
                    conclusion.append(hrDeficient).append(enter);
                }
            } else if (purity < 0.20) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.LOW_PURITY)) {
                    String lowPurity = sentenceLowPurity(purity, templateConclusion);
                    conclusion.append(lowPurity).append(enter);
                }
            }
        }

        conclusion.append("\n");

        return conclusion;
    }

    private static String createTumorLocationSentense(@NotNull String patientPrimaryTumorLocation,
            @NotNull List<TumorLocationConclusion> tumorLocationConclusion, @NotNull String cancerSubtype) {
        String locationTumor = Strings.EMPTY;
        for (TumorLocationConclusion locationConclusion : tumorLocationConclusion) {
            if (locationConclusion.primaryTumorLocation().equals(patientPrimaryTumorLocation) && locationConclusion.cancerSubType()
                    .equals(cancerSubtype)) {
                locationTumor = locationConclusion.tumorLocationConclusion();
            }
        }
        if (locationTumor.equals(Strings.EMPTY)) {
            LOGGER.warn("No tumor location is known");
        }
        return locationTumor + " sample showing: ";
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
