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
        if (systemUser.equals("lieke") || systemUser.equals("liekeschoenmaker")) {
            userName = " (Lieke Schoenmaker)";
        } else if (systemUser.equals("korneel") || systemUser.equals("korneelduyvesteyn")) {
            userName = " (Korneel Duyvesteyn)";
        } else if (systemUser.equals("sandra") || systemUser.equals("sandravandenbroek")) {
            userName = " (Sandra van den Broek)";
        } else if (systemUser.equals("root")) {
            userName = Strings.EMPTY;
        } else {
            userName = systemUser;
        }
        return "trained IT employee" + userName + " and trained clinical molecular biologist (KMBP)";
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
