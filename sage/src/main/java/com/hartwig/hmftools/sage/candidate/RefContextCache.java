package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.max;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.select.PanelSelector;

public class RefContextCache
{
    private final SageConfig mConfig;
    private final EvictingArray mEvictingArray;
    private final PanelSelector mPanelSelector;
    private final List<AltContext> mSavedCandidates;
    private final HotspotSelector mHotspotSelector;

    public RefContextCache(final SageConfig config, final List<SimpleVariant> hotspots, final List<BaseRegion> panel)
    {
        mConfig = config;
        mPanelSelector = new PanelSelector(panel);
        mSavedCandidates = Lists.newArrayList();

        mHotspotSelector = new HotspotSelector(hotspots);

        final Consumer<RefContext> evictionHandler = (refContext) -> processAltContexts(refContext);

        int minDefaultReadLength = max(config.getReadLength(), DEFAULT_READ_LENGTH);
        mEvictingArray = new EvictingArray(minDefaultReadLength, evictionHandler);
    }

    public PanelSelector panelSelector() { return mPanelSelector; }

    public RefContext getOrCreateRefContext(final String chromosome, int position)
    {
        return mEvictingArray.getOrCreateRefContext(position, aLong -> new RefContext(chromosome, position));
    }

    public List<AltContext> altContexts()
    {
        mEvictingArray.evictAll();
        Collections.sort(mSavedCandidates);
        return mSavedCandidates;
    }

    private void processAltContexts(final RefContext refContext)
    {
        Collection<AltContext> altContexts = refContext.altContexts();

        if(altContexts == null)
            return;

        for(AltContext altContext : altContexts)
        {
            if(!hasValidDnaBases(altContext))
                continue;

            if(!passesTumorHardLimits(altContext))
                continue;

            altContext.selectCandidates();

            if(altContext.hasValidCandidate())
                mSavedCandidates.add(altContext);

            if(altContext.hasSecondCandidate())
                mSavedCandidates.add(altContext.secondCandidate());
        }
    }

    private boolean hasValidDnaBases(final AltContext altContext)
    {
        for(int i = 0; i < altContext.ref().length(); i++)
        {
            if(!Nucleotides.isValidDnaBase(altContext.ref().charAt(i)))
                return false;
        }

        return true;
    }

    public boolean passesTumorHardLimits(final AltContext altContext)
    {
        if(mHotspotSelector.isHotspot(altContext))
            return true;

        return altContext.rawAltSupport() >= mConfig.Filter.HardMinTumorRawAltSupport;
    }
}
