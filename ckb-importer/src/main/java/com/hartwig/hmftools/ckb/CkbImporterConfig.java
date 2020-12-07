package com.hartwig.hmftools.ckb;

import java.io.File;
import java.nio.file.Files;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CkbImporterConfig {

    String CLINICAL_TRIALS_DIR = "clinical_trials_dir";
    String DRUG_CLASSES_DIR = "drug_classes_dir";
    String DRUGS_DIR = "drugs_dir";
    String GENES_DIR = "genes_dir";
    String GLOBAL_THERAPY_APPROVAL_STATUSES_DIR = "global_therapy_approval_statuses_dir";
    String INDICATIONS_DIR = "indications_dir";
    String MOLECULAR_PROFILES_DIR = "molecular_profiles_dir";
    String REFERENCES_DIR = "references_dir";
    String THERAPIES_DIR = "therapies_dir";
    String TREATMENT_APPROACHES_DIR = "treatment_approaches_dir";
    String VARIANTS_DIR = "variants_dir";

    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        //TODO: implement as 1 dir and not as seperate dirs

        options.addOption(CLINICAL_TRIALS_DIR, true, "Path towards the directory holding the clincial trials data");
        options.addOption(DRUG_CLASSES_DIR, true, "Path towards the directory holding the drugs classes data");
        options.addOption(DRUGS_DIR, true, "Path towards the directory holding the drugs data");
        options.addOption(GENES_DIR, true, "Path towards the directory holding the genes data");
        options.addOption(GLOBAL_THERAPY_APPROVAL_STATUSES_DIR,
                true,
                "Path towards the directory holding the LIMS global therapy approval statuses");
        options.addOption(INDICATIONS_DIR, true, "Path towards the directory holding the indications data");
        options.addOption(MOLECULAR_PROFILES_DIR, true, "Path towards the directory holding the molecular profiles data");
        options.addOption(REFERENCES_DIR, true, "Path towards the directory holding the references data");
        options.addOption(THERAPIES_DIR, true, "Path towards the directory holding the therapies data");
        options.addOption(TREATMENT_APPROACHES_DIR, true, "Path towards the directory holding the treatment approaches data");
        options.addOption(VARIANTS_DIR, true, "Path towards the directory holding the variants data");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");

        return options;
    }

    @NotNull
    String clinicalTrialsDir();

    @NotNull
    String drugsClassesDir();

    @NotNull
    String drugsDir();

    @NotNull
    String genesDir();

    @NotNull
    String globalTherpyApprovalStatusesDir();

    @NotNull
    String IndicationsDir();

    @NotNull
    String molecularProfilesDir();

    @NotNull
    String referencesDir();

    @NotNull
    String therapiesDir();

    @NotNull
    String treatmentApproachesDir();

    @NotNull
    String variantsDir();

    @NotNull
    static CkbImporterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        return ImmutableCkbImporterConfig.builder()
                .clinicalTrialsDir(nonOptionalDir(cmd, CLINICAL_TRIALS_DIR))
                .drugsClassesDir(nonOptionalDir(cmd, DRUG_CLASSES_DIR))
                .drugsDir(nonOptionalDir(cmd, DRUGS_DIR))
                .genesDir(nonOptionalDir(cmd, GENES_DIR))
                .globalTherpyApprovalStatusesDir(nonOptionalDir(cmd, GLOBAL_THERAPY_APPROVAL_STATUSES_DIR))
                .IndicationsDir(nonOptionalDir(cmd, INDICATIONS_DIR))
                .molecularProfilesDir(nonOptionalDir(cmd, MOLECULAR_PROFILES_DIR))
                .referencesDir(nonOptionalDir(cmd, REFERENCES_DIR))
                .therapiesDir(nonOptionalDir(cmd, THERAPIES_DIR))
                .treatmentApproachesDir(nonOptionalDir(cmd, TREATMENT_APPROACHES_DIR))
                .variantsDir(nonOptionalDir(cmd, VARIANTS_DIR))
                .build();
    }

    @NotNull
    static String nonOptionalDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value) || !pathIsDirectory(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing directory: " + value);
        }

        return value;
    }

    @NotNull
    static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new ParseException("Parameter must be provided: " + param);
        }

        return value;
    }

    static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }
}
