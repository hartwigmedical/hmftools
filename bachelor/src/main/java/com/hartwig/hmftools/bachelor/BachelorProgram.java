package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.HmfGenomeRegion;

public class BachelorProgram {
    private final String mName;
    private final Predicate<VariantModel> mVcfProcessor;
    private final Predicate<VariantModel> mWhiteList;
    private final Predicate<GeneCopyNumber> mCopyNumberProcessor;
    private final Predicate<HmfGenomeRegion> mDisruptionProcessor;

    private final List<String> RequiredEffects;
    private final List<String> PanelTranscripts;

    // add white and blacklist criteria for this program

    BachelorProgram(final String name, final Predicate<VariantModel> vcfProcessor, Predicate<VariantModel> whitelist,
            final Predicate<GeneCopyNumber> copyNumberProcessor, final Predicate<HmfGenomeRegion> disruptionProcessor,
            final List<String> requiredEffects, final List<String> panelTranscripts) {

        mName = name;
        mVcfProcessor = vcfProcessor;
        mWhiteList = whitelist;
        mCopyNumberProcessor = copyNumberProcessor;
        mDisruptionProcessor = disruptionProcessor;
        this.RequiredEffects = requiredEffects;
        this.PanelTranscripts = panelTranscripts;
    }

    public String name() {
        return mName;
    }

    public Predicate<VariantModel> vcfProcessor() { return mVcfProcessor; }
    public Predicate<VariantModel> whitelist() { return mWhiteList; }

    public Predicate<GeneCopyNumber> copyNumberProcessor() {
        return mCopyNumberProcessor;
    }

    public Predicate<HmfGenomeRegion> disruptionProcessor() {
        return mDisruptionProcessor;
    }

    public List<String> requiredEffects() {
        return RequiredEffects;
    }

    public List<String> panelTranscripts() {
        return PanelTranscripts;
    }
}
