package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ExtractForNumbers {

    private ExtractForNumbers() {

    }

    @NotNull
    public static String extractAllForNumbers(boolean isQCFail, double purity, boolean hasRealiablePurity, QCFailReason reason) {
        String formNumber = Strings.EMPTY;
        if (!isQCFail) {
            formNumber = succesForNumber(purity, hasRealiablePurity);
        } else {
            formNumber = extractQCFailForNumbers(reason);
        }
        return formNumber;
    }

    @NotNull
    private static String succesForNumber(double purity, boolean hasRealiablePurity) {
        String formNumber = Strings.EMPTY;

        if (purity < 0.20 || !hasRealiablePurity) {
            formNumber = ForNumber.FOR_103.display();
        } else {
            formNumber = ForNumber.FOR_080.display();

        }
        return formNumber;
    }

    @NotNull
    public static String extractSuccesfullForNumbers(double purity, boolean hasRealiablePurity) {
        return succesForNumber(purity, hasRealiablePurity);
    }

    @NotNull
    public static String extractQCFailForNumbers(QCFailReason reason) {
        return reason.usingForNumber();
    }

}
