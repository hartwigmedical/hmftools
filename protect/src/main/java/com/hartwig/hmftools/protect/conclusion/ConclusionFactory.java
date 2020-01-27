package com.hartwig.hmftools.protect.conclusion;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.common.ReportableGainLoss;
import com.hartwig.hmftools.protect.common.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.common.ReportableVariantAnalysis;
import com.hartwig.hmftools.protect.report.chord.ChordStatus;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionFactory {
    private static final Logger LOGGER = LogManager.getLogger(ConclusionFactory.class);

    private ConclusionFactory() {

    }

    public static StringBuilder createConclusion(@NotNull String patientPrimaryTumorLocation, int tumorMTL, double tumorMTB,
            double tumorMSI, double chordScore, @NotNull List<ReportableGeneFusion> geneFusions,
            @NotNull List<ReportableGainLoss> geneCopyNumbers, @NotNull ReportableVariantAnalysis reportableVariantAnalysis,
            @NotNull List<TemplateConclusion> templateConclusionList, double purity,
            @NotNull List<TumorLocationConclusion> tumorLocationConclusion, @NotNull String cancerSubtype,
            @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull Map<String, TemplateConclusion> MapTemplateConclusion) {

        StringBuilder conclusion = new StringBuilder();
        String enter = " <enter> ";
        String startRow = "- ";

        String textTumorLocation = createTumorLocationSentense(patientPrimaryTumorLocation, tumorLocationConclusion, cancerSubtype);
        conclusion.append(textTumorLocation).append(enter);

        if (ChordStatus.formChord(chordScore).equals(ChordStatus.HR_DEFICIENT)) {
            String HRD = ChordStatus.HR_DEFICIENT.toString().replace("_", "-").toLowerCase();
            for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                String keyTemplate = entry.getKey().toLowerCase();
                TemplateConclusion templateConclusion = entry.getValue();
                if (HRD.equals(keyTemplate.toLowerCase())) {
                    String hrDeficient = sentenseHrDeficient(chordScore, templateConclusion);
                    conclusion.append(startRow).append(hrDeficient).append(enter);
                }
            }
        } else if (TumorMutationalStatus.fromLoad(tumorMTL).equals(TumorMutationalStatus.HIGH)) {
            //  String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
            //  conclusion.append(startRow).append(highMTL).append(enter);
        } else if (MicrosatelliteStatus.fromIndelsPerMb(tumorMSI).equals(MicrosatelliteStatus.MSI)) {
            // String highMSI = sentenseHighMSI(tumorMSI, templateConclusion);
            //  conclusion.append(startRow).append(highMSI).append(enter);

        } else if (purity < 0.20) {
            //  String lowPurity = sentenceLowPurity(purity, templateConclusion);
            //  conclusion.append(startRow).append(lowPurity).append(enter);
        } else if (geneFusions.size() >= 1) {
            for (ReportableGeneFusion fusion : geneFusions) {
                String startFusion = fusion.geneStart().toLowerCase() + " fusion";
                String endFusion = fusion.geneEnd().toLowerCase() + " fusion";
                for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                    String keyTemplate = entry.getKey().toLowerCase();
                    TemplateConclusion templateConclusion = entry.getValue();
                    if (keyTemplate.equals(startFusion)) {
                        String fusionConclusion = sentenseFusion(fusion, templateConclusion);
                        conclusion.append(startRow).append(fusionConclusion).append(enter);
                    } else if (keyTemplate.equals(endFusion)) {
                        sentenseFusion(fusion, templateConclusion);
                        String fusionConclusion = sentenseFusion(fusion, templateConclusion);
                        conclusion.append(startRow).append(fusionConclusion).append(enter);
                    }
                }
            }
        } else if (geneCopyNumbers.size() >= 1) {
            //TODO: use only reportable geneCopyNumbers instead of all copyNumbers
            for (ReportableGainLoss gainLoss : geneCopyNumbers) {
                String interpetation = Strings.EMPTY;
                if (gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS) {
                    interpetation = "loss";
                } else if (gainLoss.interpretation() == CopyNumberInterpretation.GAIN) {
                    interpetation = "amplification";
                }

                for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                    String keyTemplate = entry.getKey().toLowerCase();
                    TemplateConclusion templateConclusion = entry.getValue();

                    if (interpetation.equals("loss")) {
                        if (keyTemplate.contains(gainLoss.gene()) && keyTemplate.contains(interpetation)) {
                            String sentenceConclusion = templateConclusion.summaryTextStatement();
                            conclusion.append(startRow).append(sentenceConclusion).append(enter);
                        }
                    } else if (interpetation.equals("amplification") && keyTemplate.contains(interpetation)) {
                        String sentenceConclusion = templateConclusion.summaryTextStatement();
                        conclusion.append(startRow).append(sentenceConclusion).append(enter);

                    }
                }
            }

        } else if (reportableHomozygousDisruptions.size() >= 1) {

        } else if (conclusion.toString().endsWith("sample showing: <enter> ")) { // Must be the last if statement)
            // conclusion.append(startRow).append(templateConclusion.summaryTextStatement());

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

    private static String sentenseFusion(@NotNull ReportableGeneFusion fusion, @NotNull TemplateConclusion templateConclusion) {
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("(XXXX)", "(" + fusion.geneStart() + " - " + fusion.geneEnd() + ")");
        return sentence;

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
