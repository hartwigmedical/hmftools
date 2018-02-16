package com.hartwig.hmftools.svannotation.analysis;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import static java.lang.Math.abs;

public class SvClusteringConfig {

    private int mClusterBaseDistance;
    private String mOutputCsvPath;
    private boolean mUseCombinedOutputFile;

    public static final int DEFAULT_BASE_DISTANCE = 100000;

    public SvClusteringConfig()
    {
        mClusterBaseDistance = DEFAULT_BASE_DISTANCE;
        mOutputCsvPath = "";
        mUseCombinedOutputFile = false;
    }

    public void setBaseDistance(int distance)
    {
        if(distance == 0)
            mClusterBaseDistance = DEFAULT_BASE_DISTANCE;
        else
            mClusterBaseDistance = distance;
    }

    public void setOutputCsvPath(final String path) { mOutputCsvPath = path; }
    public final String getOutputCsvPath() { return mOutputCsvPath; }
    public int getClusterBaseDistance() { return mClusterBaseDistance; }

    public boolean getUseCombinedOutputFile() { return mUseCombinedOutputFile; }
    public void setUseCombinedOutputFile(boolean toggle) { mUseCombinedOutputFile = true; }

}
