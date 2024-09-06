package com.hartwig.hmftools.cup.cli;

import com.hartwig.hmftools.cup.prep.CuppaDataPrep;

public class PrepPlusPredictionMain {
    public static void main(String[] args) {
        try {
            CuppaDataPrep.main(args);
            PredictionRunner.main(args);
            System.exit(0);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
