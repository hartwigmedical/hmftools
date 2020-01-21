package com.hartwig.hmftools.protect.conclusion;

import org.apache.logging.log4j.util.Strings;

public class ConclusionFactory {

    private ConclusionFactory() {

    }

    public static String createConclusion() {
        String conclusion = Strings.EMPTY;
        String enter = "<enter>";

        String textTumorLocation = createTumorLocationSentense();
        conclusion = textTumorLocation + enter;
        return conclusion;
    }

    private static String createTumorLocationSentense(){
        String tumorLocation= Strings.EMPTY;;
        return tumorLocation + "sample showing";
    }
}
