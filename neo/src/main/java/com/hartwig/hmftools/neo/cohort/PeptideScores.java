package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.cohort.StatusResults.STATUS_MAX;

public class PeptideScores
{
    public final String Peptide;

    public final double[] Affinity;
    public final double[] Presentation;

    public PeptideScores(final String peptide, double initAffinity, double initPresentation)
    {
        Peptide = peptide;

        Affinity = new double[STATUS_MAX];
        Presentation = new double[STATUS_MAX];

        for(int i = 0; i < STATUS_MAX; ++i)
        {
            Affinity[i] = initAffinity;
            Presentation[i] = initPresentation;
        }
    }
}
