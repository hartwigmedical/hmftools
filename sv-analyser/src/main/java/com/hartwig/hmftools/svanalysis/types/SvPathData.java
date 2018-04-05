package com.hartwig.hmftools.svanalysis.types;

import java.util.HashMap;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;

public class SvPathData {

    private int mId;
    private String mDesc;

    private List<SvClusterData> mLinkedVariants;
    private List<Integer> mLinkedLengths;
    private int mConsistencyValue;
    private HashMap<Integer, Integer> mClusterCounts;
    private SvUtilities mUtils;

    SvPathData(int id, final String desc, final SvUtilities utils)
    {
        mId = id;
        mDesc = desc;
        mUtils = utils;
        mLinkedVariants = Lists.newArrayList();
        mConsistencyValue = 0;
        mClusterCounts = new HashMap();
    }

    public int getId() { return mId; }

    public String getDesc() { return mDesc; }
    public void setDesc(String desc) { mDesc = desc; }

    public List<SvClusterData> getLinkedVariants() { return mLinkedVariants; }
    public List<Integer> getLinkedLengths() { return mLinkedLengths; }

    void addVariant(final SvClusterData var, int linkedLength)
    {
        mLinkedVariants.add(var);
        mLinkedLengths.add(linkedLength);
    }

    void addVariant(final SvClusterData var)
    {
        final SvClusterData linkingVar = mLinkedVariants.isEmpty() ? null : mLinkedVariants.get(mLinkedVariants.size() - 1);

        mLinkedVariants.add(var);

        if(linkingVar == null)
            return;

//        int linkedLength = mUtils.getProximity()
//        mLinkedLengths.add(linkedLength);
    }


}
