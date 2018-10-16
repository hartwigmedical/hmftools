package com.hartwig.hmftools.svanalysis.types;

import java.util.List;

import com.google.common.collect.Lists;

public class SvGeneData {

    private String mSampleId;
    private String mChromosome;
    private String mGene;
    private String mImpact;
    private String mDriver;
    private String mDriverType;
    private long mStartPosition;
    private long mEndPosition;
    private long mStartCNRegion;
    private long mEndCNRegion;
    private String mStartRegionType;
    private String mEndRegionType;

    private List<SvVarData> mStartSvList;
    private List<SvVarData> mEndSvList;

    public static String DRIVER_TYPE_TSG = "TSG";
    public static String DRIVER_TYPE_ONCOGENE = "ONCO";
    public static String DRIVER_TYPE_FUSION = "FUSION";

    public static String DRIVER_DEL = "Del";

    public SvGeneData(final String[] csvData)
    {
        mSampleId = csvData[0];
        mGene = csvData[1];
        mImpact = csvData[2];
        mDriver = csvData[3];
        mDriverType = csvData[5];
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
    public final String impact() {return mImpact; }
    public final String driver() {return mDriver; }
    public final String driverType() {return mDriverType; }
    public long startPosition() {return mStartPosition; }
    public long endPosition() {return mEndPosition; }
    public long startCNRegion() {return mStartCNRegion; }
    public long endCNRegion() {return mEndCNRegion + 1; } // since CN region ends just before the SV position
    public final String startRegionType() {return mStartRegionType; }
    public final String endRegionType() {return mEndRegionType; }

    public void addSvData(final SvVarData var, boolean isStart)
    {
        if(isStart)
            mStartSvList.add(var);
        else
            mEndSvList.add(var);
    }

    public final List<SvVarData> getStartSvList() { return mStartSvList; }
    public final List<SvVarData> getEndSvList() { return mEndSvList; }

}
