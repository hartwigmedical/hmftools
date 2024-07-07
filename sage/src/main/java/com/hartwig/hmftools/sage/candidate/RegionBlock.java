package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.ceil;
import static java.lang.String.format;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.select.PanelSelector;
import com.hartwig.hmftools.sage.select.ReadPanelStatus;

public class RegionBlock extends BaseRegion
{
    private final ReadPanelStatus mPanelStatus;
    private final Set<Integer> mHotspotPositions;
    private final boolean mApplyEventPenalty;
    private final int mDepthLimit;

    private int mDepth;

    public RegionBlock(final int posStart, final int posEnd, final ReadPanelStatus panelStatus, final Set<Integer> hotspotPositions,
            final boolean applyEventPenalty, final int depthLimit)
    {
        super(posStart, posEnd);
        mPanelStatus = panelStatus;
        mHotspotPositions = hotspotPositions;
        mApplyEventPenalty = applyEventPenalty;
        mDepthLimit = depthLimit;
    }

    public ReadPanelStatus panelStatus() { return mPanelStatus; }

    public int depthLimit() { return mDepthLimit; }
    public int currentDepth() { return mDepth; }
    public void incrementDepth() { ++mDepth; }
    public boolean depthLimitReached() { return mDepth >= mDepthLimit; }

    public boolean isHotspot(int position) { return mHotspotPositions.contains(position); }
    public boolean coversHotspot(int readStart, int readEnd)
    {
        return mHotspotPositions.stream().anyMatch(x -> positionWithin(x, readStart, readEnd)); }

    public Set<Integer> hotspotPositions() { return mHotspotPositions; }

    public boolean applyEventPenalty() { return mApplyEventPenalty; }

    public String toString()
    {
        return format("%s depth(%d/%d) hotspots(%d) %s",
                super.toString(), mDepth, mDepthLimit, mHotspotPositions.size(), !mApplyEventPenalty ? "no-event-penalty" : "");
    }

    public static List<RegionBlock> buildRegionBlocks(
            int blockSize, final ChrBaseRegion region, final PanelSelector panelSelector, final List<SimpleVariant> regionHotspots,
            int maxReadDepthPanel, int maxReadDepthNonPanel)
    {
        boolean isMitochonrial = MitochondrialChromosome.contains(region.Chromosome);

        int blockCount = (int)ceil(region.length() / blockSize);
        List<RegionBlock> regionBlocks = Lists.newArrayListWithCapacity(blockCount);

        for(int blockPos = region.start(); blockPos <= region.end(); blockPos += blockSize)
        {
            int nextBlockStart = blockPos + blockSize;

            int blockStart = blockPos;
            int blockEnd = region.end() - nextBlockStart < blockSize / 2 ? region.end() : nextBlockStart - 1;

            Set<Integer> hotspotPositions = regionHotspots.stream()
                    .filter(x -> positionWithin(x.Position, blockStart, blockEnd))
                    .map(x -> x.Position).collect(Collectors.toSet());

            ReadPanelStatus panelStatus = panelSelector.panelStatus(blockStart, blockEnd);

            int depthLimit = !isMitochonrial && panelStatus == ReadPanelStatus.OUTSIDE_PANEL ? maxReadDepthNonPanel : maxReadDepthPanel;

            boolean applyEventPenalty = !HlaCommon.overlaps(region.Chromosome, blockStart, blockEnd);

            RegionBlock regionBlock = new RegionBlock(blockStart, blockEnd, panelStatus, hotspotPositions, applyEventPenalty, depthLimit);

            regionBlocks.add(regionBlock);

            if(blockEnd == region.end())
                break;
        }

        return regionBlocks;
    }
}
