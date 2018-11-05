package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import nl.hartwigmedicalfoundation.bachelor.HotspotLocation;

public class BachelorProgram
{
    private final String mName;
    private final Predicate<VariantModel> mVcfProcessor;
    private final Predicate<VariantModel> mWhiteList;

    private final List<String> mRequiredEffects;
    private final List<String> mPanelTranscripts;
    private final List<HotspotLocation> mHotspots;

    BachelorProgram(final String name, final Predicate<VariantModel> vcfProcessor, Predicate<VariantModel> whitelist,
            final List<String> requiredEffects, final List<String> panelTranscripts, final List<HotspotLocation> hotspots)
    {
        mName = name;
        mVcfProcessor = vcfProcessor;
        mWhiteList = whitelist;
        mRequiredEffects = requiredEffects;
        mPanelTranscripts = panelTranscripts;
        mHotspots = hotspots;
    }

    public String name() {
        return mName;
    }

    public Predicate<VariantModel> vcfProcessor() { return mVcfProcessor; }
    public Predicate<VariantModel> whitelist() { return mWhiteList; }

    public List<String> requiredEffects() { return mRequiredEffects; }
    public List<String> panelTranscripts() { return mPanelTranscripts; }
    public List<HotspotLocation> hotspots() { return mHotspots; }
}
