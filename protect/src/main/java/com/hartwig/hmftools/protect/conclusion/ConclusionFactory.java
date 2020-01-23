package com.hartwig.hmftools.protect.conclusion;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.common.ReportableGainLoss;
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
            double chordScore, @NotNull List<ReportableGeneFusion> geneFusions, @NotNull List<ReportableGainLoss> geneCopyNumbers,
            @NotNull List<? extends Variant> passSomaticVariants, @NotNull List<TemplateConclusion> templateConclusionList, double purity,
            @NotNull List<TumorLocationConclusion> tumorLocationConclusion, @NotNull String cancerSubtype) {

        StringBuilder conclusion = new StringBuilder();
        String enter = " <enter> ";
        String startRow = "- ";

        LOGGER.info(geneCopyNumbers);


        String textTumorLocation = createTumorLocationSentense(patientPrimaryTumorLocation, tumorLocationConclusion, cancerSubtype);
        conclusion.append(textTumorLocation).append(enter);

        for (TemplateConclusion templateConclusion : templateConclusionList) {
            if (TumorMutationalStatus.fromLoad(tumorMTL).equals(TumorMutationalStatus.HIGH)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.HIGH_MTL)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(startRow).append(highMTL).append(enter);
                }
            } else if (MicrosatelliteStatus.fromIndelsPerMb(tumorMSI).equals(MicrosatelliteStatus.MSI)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.MSI)) {
                    String highMSI = sentenseHighMSI(tumorMSI, templateConclusion);
                    conclusion.append(startRow).append(highMSI).append(enter);
                }
            } else if (ChordStatus.formChord(chordScore).equals(ChordStatus.HR_DEFICIENT)) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.HR_DEFICIENT)) {
                    String hrDeficient = sentenseHrDeficient(chordScore, templateConclusion);
                    conclusion.append(startRow).append(hrDeficient).append(enter);
                }
            } else if (purity < 0.20) {
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.LOW_PURITY)) {
                    String lowPurity = sentenceLowPurity(purity, templateConclusion);
                    conclusion.append(startRow).append(lowPurity).append(enter);
                }
            } else if (geneFusions.size() >= 1) {
                for (ReportableGeneFusion fusion: geneFusions) {
                    String geneFusionStart = fusion.geneStart() + " fusion";
                    String geneFusionEnd = fusion.geneEnd() + " fusion";
                    if (geneFusionStart.equals(templateConclusion.abberrationGeneSummary()) || geneFusionEnd.equals(templateConclusion.abberrationGeneSummary())) {
                        if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.BRAF_FUSION)){
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.ALK_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.EGFR_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.FGFR1_FUSION)){
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.FGFR2_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.FGFR3_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.FGFR4_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.MET_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.NRG1_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.NTRK2_FUSION)){
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.NTRK3_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.RET_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.ROS1_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.RSPO2_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.RSPO3_FUSION)) {
                            conclusion.append(startRow).append(templateConclusion.summaryTextStatement()).append(enter);
                        } else {
                            LOGGER.warn("No evidence text present for fusion. The fusion gene start is {} and the fusion gene end is {}.", geneFusionStart, geneFusionEnd);
                        }
                    }
                }
            }



            else if (conclusion.toString().endsWith("sample showing: <enter> ")) { // Must be the last if statement
                if (AberrationGenSummary.fromString(templateConclusion.abberrationGeneSummary()).equals(AberrationGenSummary.NONE_FOUND)) {
                    conclusion.append(startRow).append(templateConclusion.summaryTextStatement());
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
        return locationTumor + " sample showing:";
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
