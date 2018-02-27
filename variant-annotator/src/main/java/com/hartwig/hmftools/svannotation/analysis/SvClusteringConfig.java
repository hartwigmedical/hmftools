package com.hartwig.hmftools.svannotation.analysis;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import static java.lang.Math.abs;

public class SvClusteringConfig {

    private int mClusterBaseDistance;
    private String mOutputCsvFile;

    public static final int DEFAULT_BASE_DISTANCE = 1000;

    public SvClusteringConfig()
    {
        mClusterBaseDistance = DEFAULT_BASE_DISTANCE;
        mOutputCsvFile = "";
    }

    public void setBaseDistance(int distance)
    {
        if(distance == 0)
            mClusterBaseDistance = DEFAULT_BASE_DISTANCE;
        else
            mClusterBaseDistance = distance;
    }

    public void setOutputCsvFile(final String outputFile) { mOutputCsvFile = outputFile; }

    public final String getOutputCsvFile() { return mOutputCsvFile; }
    public int getClusterBaseDistance() { return mClusterBaseDistance; }

}
