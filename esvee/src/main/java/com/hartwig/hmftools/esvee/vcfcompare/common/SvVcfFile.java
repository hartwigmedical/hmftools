package com.hartwig.hmftools.esvee.vcfcompare.common;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.VariantContext;

public class SvVcfFile
{
    public final String mPath;
    public final String mLabel;

    private final SvCallerType mSvCallerType;
    private final VcfType mSourceVcfType;

    private Map<String, List<VariantBreakend>> mChrBreakendMap;
    private List<VariantBreakend> mBreakendListCache;

    public SvVcfFile(String path, String label)
    {
        mPath = path;
        mLabel = label;

        mChrBreakendMap = null;
        mBreakendListCache = null;

        mSvCallerType = SvCallerType.fromVcfPath(mPath);
        mSourceVcfType = VcfType.fromVcfPath(mPath);
    }

    public SvVcfFile loadVariants()
    {
        SV_LOGGER.info("loading {} from file({})", mLabel, mPath);

        VcfFileReader reader = new VcfFileReader(mPath);

        int variantCount = 0;
        Map<String, List<VariantBreakend>> chrBreakendMap = new HashMap<>();

        for(VariantContext variantContext : reader.iterator())
        {
            VariantBreakend variantBreakend = new VariantBreakend(variantContext, mSvCallerType, mSourceVcfType);

            String chromosome = variantContext.getContig();
            if(!chrBreakendMap.containsKey(chromosome))
            {
                chrBreakendMap.put(chromosome, new ArrayList<>());
            }
            chrBreakendMap.get(chromosome).add(variantBreakend);

            variantCount++;
        }

        mChrBreakendMap = chrBreakendMap;

        SV_LOGGER.debug("loaded {} SVs", variantCount);

        return this;
    }

    private void checkVariantsLoaded()
    {
        if(mChrBreakendMap == null)
            throw new IllegalStateException("loadVariants() not yet called");
    }

    public Map<String,List<VariantBreakend>> getVariantsAsMap()
    {
        checkVariantsLoaded();
        return mChrBreakendMap;
    }

    public List<VariantBreakend> getVariantsAsList()
    {
        checkVariantsLoaded();
        if(mBreakendListCache == null)
        {
            mBreakendListCache = new ArrayList<>();
            for(String chromosome : mChrBreakendMap.keySet())
            {
                mBreakendListCache.addAll(mChrBreakendMap.get(chromosome));
            }
        }

        return mBreakendListCache;
    }
}
