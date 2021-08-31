package com.hartwig.hmftools.patientreporter;

import java.util.Optional;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public interface PatientReport {

    @NotNull
    SampleReport sampleReport();

    @NotNull
    default String user() {
        String systemUser = System.getProperty("user.name");
        String userName = Strings.EMPTY;
        String trainedEmployee = " (trained IT employee)";
        String combinedUserName = Strings.EMPTY;
        if (systemUser.equals("lieke") || systemUser.equals("liekeschoenmaker")) {
            userName = "Lieke Schoenmaker";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("korneel") || systemUser.equals("korneelduyvesteyn")) {
            userName = "Korneel Duyvesteyn";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("sandra") || systemUser.equals("sandravandenbroek")) {
            userName = "Sandra van den Broek";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("root")) {
            combinedUserName = "automatically";
        } else {
            userName = systemUser;
            combinedUserName = userName + trainedEmployee;
        }

        if (!combinedUserName.endsWith(trainedEmployee)) {
            combinedUserName = "by " + combinedUserName;
        }

        return combinedUserName + " and checked by a trained Clinical Molecular Biologist in Pathology (KMBP)";

    }

    @NotNull
    String qsFormNumber();

    @NotNull
    Optional<String> comments();

    boolean isCorrectedReport();

    @NotNull
    String signaturePath();

    @NotNull
    String logoRVAPath();

    @NotNull
    String logoCompanyPath();
}
