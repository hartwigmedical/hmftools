package com.hartwig.hmftools.protect.conclusion;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.common.DriverInterpretation;
import com.hartwig.hmftools.protect.common.ReportableGainLoss;
import com.hartwig.hmftools.protect.common.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.common.ReportableVariant;
import com.hartwig.hmftools.protect.common.ReportableVariantAnalysis;

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
            @NotNull List<ReportableGainLoss> geneCopyNumbers, @NotNull ReportableVariantAnalysis reportableVariantAnalysis, double purity,
            @NotNull String cancerSubtype, @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull Map<String, TemplateConclusion> MapTemplateConclusion) {

        StringBuilder conclusion = new StringBuilder();
        String enter = " <enter> ";
        String startRow = "- ";

        String textTumorLocation = createTumorLocationSentense(patientPrimaryTumorLocation, cancerSubtype);
        conclusion.append(textTumorLocation).append(enter);

        if (chordScore > 0.5) {
            String HRD = "HR-deficient";
            for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                TemplateConclusion templateConclusion = entry.getValue();
                if (HRD.equals("HR-deficient")) {
                    String hrDeficient = sentenceHrDeficient(chordScore, templateConclusion);
                    conclusion.append(startRow).append(hrDeficient).append(enter);
                }
            }
        }

        String melanoma = cancerSubtype.toLowerCase();
        String lung = patientPrimaryTumorLocation.toLowerCase();

        if (TumorMutationalStatus.fromLoad(tumorMTL).equals(TumorMutationalStatus.HIGH)) {
            TumorMutationalStatus TML = TumorMutationalStatus.fromLoad(tumorMTL);
            String TMLValue = Strings.EMPTY;
            if (TML.equals(TumorMutationalStatus.HIGH)) {
                TMLValue = "high tmb";
            } else {
                TMLValue = "low tmb";
            }
            for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                String keyTemplate = entry.getKey().toLowerCase();
                TemplateConclusion templateConclusion = entry.getValue();
                if (TMLValue.contains(keyTemplate) && keyTemplate.contains(melanoma)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(startRow).append(highMTL).append(enter);
                } else if (TMLValue.contains(keyTemplate) && keyTemplate.contains(lung)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(startRow).append(highMTL).append(enter);
                } else if (TMLValue.equals(keyTemplate) && !keyTemplate.contains(melanoma)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(startRow).append(highMTL).append(enter);
                } else if (TMLValue.equals(keyTemplate) && !keyTemplate.contains(lung)) {
                    String highMTL = sentenseHighMTL(tumorMTL, tumorMTB, templateConclusion);
                    conclusion.append(startRow).append(highMTL).append(enter);
                }
            }
        }

        if (MicrosatelliteStatus.fromIndelsPerMb(tumorMSI).equals(MicrosatelliteStatus.MSI)) {
            MicrosatelliteStatus MSI = MicrosatelliteStatus.fromIndelsPerMb(tumorMSI);
            String MSIValue = Strings.EMPTY;
            if (MSI.equals(MicrosatelliteStatus.MSI)) {
                MSIValue = "msi";
            } else {
                MSIValue = "mss";
            }
            for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                String keyTemplate = entry.getKey().toLowerCase();
                TemplateConclusion templateConclusion = entry.getValue();
                if (MSIValue.equals(keyTemplate.toLowerCase())) {
                    String highMSI = sentenseHighMSI(tumorMSI, templateConclusion);
                    conclusion.append(startRow).append(highMSI).append(enter);
                }
            }
        }

        if (purity < 0.20) {
            String purityValue = "LOW PURITY".toLowerCase();
            for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                String keyTemplate = entry.getKey().toLowerCase();
                TemplateConclusion templateConclusion = entry.getValue();
                if (purityValue.equals(keyTemplate.toLowerCase())) {
                    String lowPurity = sentenceLowPurity(purity, templateConclusion);
                    conclusion.append(startRow).append(lowPurity).append(enter);
                }
            }
        }

        if (geneFusions.size() >= 1) {
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
        }

        if (geneCopyNumbers.size() >= 1) {
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

                    if (keyTemplate.contains(gainLoss.gene().toLowerCase()) && keyTemplate.contains(interpetation)) {
                        String sentenceConclusion = templateConclusion.summaryTextStatement();
                        conclusion.append(startRow).append(sentenceConclusion).append(enter);
                    }
                }
            }
        }

        if (reportableHomozygousDisruptions.size() >= 1) {
            for (ReportableHomozygousDisruption reportableHomozygousDisruption : reportableHomozygousDisruptions) {
                for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                    String keyTemplate = entry.getKey().toLowerCase();
                    TemplateConclusion templateConclusion = entry.getValue();

                    if (keyTemplate.contains(reportableHomozygousDisruption.gene().toLowerCase()) && keyTemplate.contains("inactivation")) {
                        String sentenceConclusion = templateConclusion.summaryTextStatement();
                        conclusion.append(startRow).append(sentenceConclusion).append(enter);
                    }
                }
            }
        }

        if (reportableVariantAnalysis.variantsToReport().size() >= 1) {
            List<ReportableVariant> reportableVariants = reportableVariantAnalysis.variantsToReport();
            for (ReportableVariant variant : reportableVariants) {
                DriverInterpretation interpretation = DriverInterpretation.interpret(variant.driverLikelihood());

                if (interpretation == DriverInterpretation.HIGH) {
                    for (Map.Entry<String, TemplateConclusion> entry : MapTemplateConclusion.entrySet()) {
                        String keyTemplate = entry.getKey().toLowerCase();
                        TemplateConclusion templateConclusion = entry.getValue();
                        if (keyTemplate.contains(variant.gene().toLowerCase()) && keyTemplate.contains("mutation")) {
                            String sentenceConclusion = templateConclusion.summaryTextStatement();
                            conclusion.append(startRow).append(sentenceConclusion).append(enter);
                        }
                    }
                }
            }
        }

        if (conclusion.toString().endsWith("sample showing: <enter> ")) { // Must be the last if statement)

            //   conclusion.append(startRow).append(templateConclusion.summaryTextStatement());
        }

        conclusion.append("\n");

        return conclusion;
    }

    private static String createTumorLocationSentense(@NotNull String patientPrimaryTumorLocation, @NotNull String cancerSubtype) {
        String tumorLocation = Strings.EMPTY;
        if (patientPrimaryTumorLocation.equals(Strings.EMPTY)) {
            LOGGER.warn("No tumor location is known of patient!");
            tumorLocation = "undetermined";
        } else if (!patientPrimaryTumorLocation.equals(Strings.EMPTY) && !cancerSubtype.equals(Strings.EMPTY)) {
            tumorLocation = patientPrimaryTumorLocation + " (" + cancerSubtype + ")";
        } else if (!patientPrimaryTumorLocation.equals(Strings.EMPTY) && cancerSubtype.equals(Strings.EMPTY)) {
            tumorLocation = patientPrimaryTumorLocation;
        }
        return tumorLocation + " cancer sample showing:";
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

    private static String sentenceHrDeficient(double chordScore, @NotNull TemplateConclusion templateConclusion) {
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("a high CHORD score (XXX) ", "a high CHORD score (" + String.format("%.2f", chordScore) + ") ");
        return sentence;
    }

    private static String sentenceLowPurity(double purity, @NotNull TemplateConclusion templateConclusion) {
        String purityPercentage = new DecimalFormat("#'%'").format(purity * 100);
        String sentence = templateConclusion.summaryTextStatement();
        sentence = sentence.replace("Due to the low tumor purity (XXX) ", "Due to the low tumor purity (" + purityPercentage + ") ");
        return sentence;
    }
}
