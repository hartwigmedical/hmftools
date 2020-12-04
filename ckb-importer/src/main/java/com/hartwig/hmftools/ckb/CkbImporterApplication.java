package com.hartwig.hmftools.ckb;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrial;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbImporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(CkbImporterApplication.class);
    public static final String VERSION = CkbImporterApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running CKB importer v{}", VERSION);

        Options options = CkbImporterConfig.createOptions();

        CkbImporterConfig config = null;
        try {
            config = CkbImporterConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("CKB Importer", options);
            System.exit(1);
        }

        readJsonData(config);
        LOGGER.info("Complete!");
    }

    private static void readJsonData(@NotNull CkbImporterConfig config) throws IOException {
        List<ClinicalTrial> clinicalTrials = readingClinicalTrial(config.clinicalTrialsDir());
        readingDrugsClasses(config.drugsClassesDir());
        readingDrugs(config.drugsDir());
        readingGenes(config.genesDir());
        readingGlobalTherapyApprovalStatuses(config.globalTherpyApprovalStatusesDir());
        readingIndications(config.IndicationsDir());
        readingMolecularProfiles(config.molecularProfilesDir());
        readingReferences(config.referencesDir());
        readingTherapies(config.therapiesDir());
        readingTreatmentApproaches(config.treatmentApproachesDir());
        readingVariants(config.variantsDir());

        //implement in 1 object
    }

    @NotNull
    private static List<ClinicalTrial> readingClinicalTrial(@NotNull String clinicalTrialDir) throws IOException {
        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        File[] filesClinicalTrials = new File(clinicalTrialDir).listFiles();
        if (filesClinicalTrials != null) {
            for (File clinicalTrial : filesClinicalTrials) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(clinicalTrial));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject clinicalTrialsEntryObject = parser.parse(reader).getAsJsonObject();
                    //clinicalTrials.add();
                }
            }
        }
        return clinicalTrials;
    }

    private static void readingDrugsClasses(@NotNull String drugsClassesDir) throws IOException {

    }

    private static void readingDrugs(@NotNull String drugsDir) throws IOException {

    }

    private static void readingGenes(@NotNull String genesDir) throws IOException {

    }

    private static void readingGlobalTherapyApprovalStatuses(@NotNull String globalTherpyApprovalStatusesDir) throws IOException {

    }

    private static void readingIndications(@NotNull String indicationsDir) throws IOException {

    }

    private static void readingMolecularProfiles(@NotNull String molecularProfilesDir) throws IOException {

    }

    private static void readingReferences(@NotNull String referencesDir) throws IOException {

    }

    private static void readingTreatmentApproaches(@NotNull String therapiesDir) throws IOException {

    }

    private static void readingTherapies(@NotNull String treatmentApproachesDir) throws IOException {

    }

    private static void readingVariants(@NotNull String variantsDir) throws IOException {

    }
}
