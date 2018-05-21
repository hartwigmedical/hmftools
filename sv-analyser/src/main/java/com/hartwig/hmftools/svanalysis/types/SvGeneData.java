package com.hartwig.hmftools.svanalysis.types;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.analysis.SvCluster;

public class SvGeneData {

    private String mSampleId;
    private String mChromosome;
    private String mGene;
    private String mDriveType;
    private long mStartPosition;
    private long mEndPosition;
    private long mStartCNRegion;
    private long mEndCNRegion;
    private String mStartRegionType;
    private String mEndRegionType;

    private List<SvClusterData> mStartSvList;
    private List<SvClusterData> mEndSvList;

    public static String DRIVER_TYPE_TSG = "TSG";
    public static String DRIVER_TYPE_ONCOGENE = "ONCO";
    public static String DRIVER_TYPE_FUSION = "FUSION";

    public SvGeneData(final String[] csvData)
    {
        mSampleId = csvData[0];
        mGene = csvData[1];
        mDriveType = csvData[5];
        mChromosome = csvData[10];
        mStartPosition = Long.parseLong(csvData[11]);
        mEndPosition = Long.parseLong(csvData[12]);
        mStartCNRegion = Long.parseLong(csvData[16]);
        mEndCNRegion = Long.parseLong(csvData[17]);
        mStartRegionType = csvData[18];
        mEndRegionType = csvData[19];

        mStartSvList = Lists.newArrayList();
        mEndSvList = Lists.newArrayList();
    }

    public final String sampleId() {return mSampleId; }
    public final String chromosome() {return mChromosome ; }
    public final String gene() {return mGene; }
    public final String driveType() {return mDriveType; }
    public long startPosition() {return mStartPosition; }
    public long endPosition() {return mEndPosition + 1; } // since CN region ends just before the SV position
    public long startCNRegion() {return mStartCNRegion; }
    public long endCNRegion() {return endCNRegion(); }
    public final String startRegionType() {return mStartRegionType; }
    public final String endRegionType() {return mEndRegionType; }

    public void addSvData(final SvClusterData var, boolean isStart)
    {
        if(isStart)
            mStartSvList.add(var);
        else
            mEndSvList.add(var);
    }

    public final List<SvClusterData> getStartSvList() { return mStartSvList; }
    public final List<SvClusterData> getEndSvList() { return mEndSvList; }

}
