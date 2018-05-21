package com.hartwig.hmftools.svanalysis.analysis;

public class SvClusteringConfig {

    private int mClusterBaseDistance;
    private String mOutputCsvPath;
    private String mSvPONFile;
    private String mFragileSiteFile;
    private String mLineElementFile;
    private String mExternalAnnotationsFile;
    private String mGeneDataFile;
    private boolean mUseCombinedOutputFile;

    public static final int DEFAULT_BASE_DISTANCE = 100000;

    public SvClusteringConfig()
    {
        mClusterBaseDistance = DEFAULT_BASE_DISTANCE;
        mFragileSiteFile = "";
        mLineElementFile = "";
        mExternalAnnotationsFile = "";
        mGeneDataFile = "";

        mOutputCsvPath = "";
        mUseCombinedOutputFile = false;
        mSvPONFile = "";
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

    public void setSvPONFile(final String filename) { mSvPONFile = filename; }
    public String getSvPONFile() { return mSvPONFile; }

    public String getFragileSiteFile() { return mFragileSiteFile;}
    public void setFragileSiteFile(final String filename) { mFragileSiteFile = filename;}

    public void setLineElementFile(final String filename) { mLineElementFile = filename; }
    public String getLineElementFile() { return mLineElementFile; }

    public void setExternalAnnotationsFile(final String filename) { mExternalAnnotationsFile = filename; }
    public String getExternalAnnotationsFile() { return mExternalAnnotationsFile; }

    public void setGeneDataFile(final String filename) { mGeneDataFile = filename; }
    public String getGeneDataFile() { return mGeneDataFile; }

    public int getClusterBaseDistance() { return mClusterBaseDistance; }

    public boolean getUseCombinedOutputFile() { return mUseCombinedOutputFile; }
    public void setUseCombinedOutputFile(boolean toggle) { mUseCombinedOutputFile = toggle; }

}
