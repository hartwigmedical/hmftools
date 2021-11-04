package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
        } else if (systemUser.equals("sandra") || systemUser.equals("sandravandenbroek") || systemUser.equals("sandravdbroek")
                || systemUser.equals("s_vandenbroek")) {
            userName = "Sandra van den Broek";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("daphne") || systemUser.equals("d_vanbeek") || systemUser.equals("daphnevanbeek")) {
            userName = "Daphne van Beek";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("root")) {
            combinedUserName = "automatically";
        } else {
            userName = systemUser;
            combinedUserName = userName + trainedEmployee;
        }

        if (combinedUserName.endsWith(trainedEmployee)) {
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

    @NotNull
    String udiDi();

    @NotNull
    List<PeachGenotype> peachGenotypes();
}
