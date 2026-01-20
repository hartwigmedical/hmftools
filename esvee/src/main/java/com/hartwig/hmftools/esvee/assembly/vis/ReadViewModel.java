package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.SvgRender.BOX_PADDING;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.INDEL_CORRECTION;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.READ_HEIGHT_PX;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static j2html.TagCreator.div;
import static j2html.TagCreator.rawHtml;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.List;
import java.util.NavigableMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffInfo;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffType;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisualiser.SegmentViewModel;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import j2html.tags.DomContent;

public final class ReadViewModel
{
    private final SupportRead mSupportRead;
    private final List<SegmentViewModel> mRefViewModel;
    private final BaseSeqViewModel mReadViewModel;

    private ReadViewModel(final SupportRead supportRead, final List<SegmentViewModel> refViewModel, final BaseSeqViewModel readViewModel)
    {
        mSupportRead = supportRead;
        mRefViewModel = refViewModel;
        mReadViewModel = readViewModel;
    }

    private static NavigableMap<Integer, List<CigarOperator>> mapCigarOperators(final List<CigarElement> cigarEls)
    {
        NavigableMap<Integer, List<CigarOperator>> mappedOps = Maps.newTreeMap();
        int baseIdx = 0;
        for(CigarElement cigarEl : cigarEls)
        {
            CigarOperator op = cigarEl.getOperator();
            int length = cigarEl.getLength();
            for(int i = 0; i < length; i++)
            {
                if(op.consumesReadBases())
                {
                    mappedOps.computeIfAbsent(baseIdx, k -> Lists.newArrayList()).add(op);
                    baseIdx++;
                }
                else
                {
                    mappedOps.computeIfAbsent(baseIdx - 1, k -> Lists.newArrayList()).add(op);
                }
            }
        }

        return mappedOps;
    }

    private static List<CigarElement> collapseMappedCigarOperators(final NavigableMap<Integer, List<CigarOperator>> mappedOps)
    {
        List<CigarElement> cigarEls = Lists.newArrayList();
        CigarOperator currentOp = null;
        int currentLength = 0;
        for(List<CigarOperator> ops : mappedOps.values())
        {
            for(CigarOperator op : ops)
            {
                if(currentOp == null)
                {
                    currentOp = op;
                    currentLength = 1;
                    continue;
                }

                if(currentOp == op)
                {
                    currentLength += 1;
                    continue;
                }

                cigarEls.add(new CigarElement(currentLength, currentOp));
                currentOp = op;
                currentLength = 1;
            }
        }

        if(currentLength > 0)
            cigarEls.add(new CigarElement(currentLength, currentOp));

        return cigarEls;
    }

    private record CorrectIndelsResult(List<CigarElement> cigarEls, int indexOffset) {}

    private static CorrectIndelsResult correctIndels(final List<SegmentViewModel> refViewModel, final SupportRead read, final List<CigarElement> cigarEls)
    {
        BaseRegion alignment = new BaseRegion(read.alignmentStart(), read.alignmentEnd());
        Boolean isBuiltForward = null;
        SegmentViewModel firstRefViewModel = refViewModel.get(0);
        SegmentViewModel lastRefViewModel = refViewModel.get(refViewModel.size() - 1);
        if(read.chromosome().equals(firstRefViewModel.chromosome()) && alignment.overlaps(firstRefViewModel.refRegion()))
            isBuiltForward = true;
        else if(read.chromosome().equals(lastRefViewModel.chromosome()) && alignment.overlaps(lastRefViewModel.refRegion()))
            isBuiltForward = false;

        if(isBuiltForward == null)
        {
            SV_LOGGER.error("Cannot find matching junction for read: {}", read.id());
            System.exit(0);
        }

        List<SequenceDiffInfo> mismatches = read.cachedRead().mismatches();
        if(mismatches == null || mismatches.isEmpty())
            return new CorrectIndelsResult(cigarEls, 0);

        NavigableMap<Integer, List<CigarOperator>> mappedCigarOps = mapCigarOperators(cigarEls);
        int indelOffset = 0;
        for(SequenceDiffInfo mismatch : mismatches)
        {
            if(mismatch.Type == SequenceDiffType.BASE)
                continue;

            if(mismatch.Type == SequenceDiffType.INSERT)
            {
                int insertLength = abs(mismatch.IndelLength);
                if(!isBuiltForward)
                    indelOffset += insertLength;

                for(int i = 0; i < insertLength; i++)
                {
                    int baseIdx = mismatch.ReadIndex + i * (isBuiltForward ? 1 : -1);
                    mappedCigarOps.get(baseIdx).set(0, I);
                }

                continue;
            }

            if(mismatch.Type == SequenceDiffType.DELETE)
            {
                int delLength = abs(mismatch.IndelLength);
                if(!isBuiltForward)
                    indelOffset -= delLength;

                int baseIdx = mismatch.ReadIndex;
                for(int i = 0; i < delLength; i++)
                    mappedCigarOps.get(baseIdx).add(D);

                continue;
            }

            if(mismatch.Type == SequenceDiffType.REPEAT)
            {
                if(mismatch.IndelLength > 0)
                {
                    int insertLength = abs(mismatch.IndelLength);
                    if(!isBuiltForward)
                        indelOffset += insertLength;

                    for(int i = 0; i < insertLength; i++)
                    {
                        int baseIdx = mismatch.ReadIndex + (i + 1) * (isBuiltForward ? 1 : -1);
                        mappedCigarOps.get(baseIdx).set(0, I);
                    }
                }
                else
                {
                    int delLength = abs(mismatch.IndelLength);
                    if(!isBuiltForward)
                        indelOffset -= delLength;

                    int baseIdx = mismatch.ReadIndex - (isBuiltForward ? 0 : 1);
                    for(int i = 0; i < delLength; i++)
                        mappedCigarOps.get(baseIdx).add(D);
                }

                continue;
            }

            // TODO(mkcmkc): Handle different mismatch types.
            if(true)
                throw new NotImplementedException("TODO");
        }

        List<CigarElement> correctedCigarEls = collapseMappedCigarOperators(mappedCigarOps);
        return new CorrectIndelsResult(correctedCigarEls, indelOffset);
    }

    public static ReadViewModel create(final List<SegmentViewModel> refViewModel, final SupportRead read, final JunctionAssembly junctionAssembly)
    {
        boolean readNegativeStrandFlag = read.orientation() == REVERSE;
        byte[] readBases = read.cachedRead().getBases();
        byte[] readBaseQuals = read.cachedRead().getBaseQuality();
        if(read.fullAssemblyOrientation() == REVERSE)
        {
            // TODO(mkcmkc): how does this interact with mismatch info?
            if(true)
                throw new NotImplementedException("TODO");

            readBases = reverseComplementBases(readBases);

            int left = 0;
            int right = readBaseQuals.length - 1;
            while(left < right)
            {
                byte temp = readBaseQuals[left];
                readBaseQuals[left] = readBaseQuals[right];
                readBaseQuals[right] = temp;
                left += 1;
                right -= 1;
            }
        }

        List<CigarElement> cigarEls = read.cachedRead().cigarElements();
        int indelOffset = 0;
        if(junctionAssembly.isForwardJunction())
        {
            for(CigarElement cigarEl : cigarEls)
            {
                if(cigarEl.getOperator() == I)
                    indelOffset += cigarEl.getLength();
                else if(cigarEl.getOperator() == D)
                    indelOffset -= cigarEl.getLength();
            }
        }

        if(INDEL_CORRECTION)
        {
            CorrectIndelsResult correctIndelsResult = correctIndels(refViewModel, read, cigarEls);
            cigarEls = correctIndelsResult.cigarEls;
            indelOffset += correctIndelsResult.indexOffset;
        }

        BaseSeqViewModel readViewModel = BaseSeqViewModel.create(
                read.fullAssemblyIndexStart() + indelOffset, cigarEls, readBases, readBaseQuals, readNegativeStrandFlag);
        readViewModel.clearSoftClips();

        return new ReadViewModel(read, refViewModel, readViewModel);
    }

    @Nullable
    public DomContent render()
    {
        BaseRegion readRegion = new BaseRegion(mReadViewModel.FirstBasePos, mReadViewModel.LastBasePos);

        int totalBoxWidth = 0;
        for(SegmentViewModel refViewModel : mRefViewModel)
            totalBoxWidth += refViewModel.viewRegion().baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        int renderCount = 0;
        for(int i = 0; i < mRefViewModel.size(); i += 1)
        {
            SegmentViewModel refEl = mRefViewModel.get(i);
            BaseRegion viewRegion = refEl.viewRegion();

            if(!viewRegion.overlaps(readRegion))
            {
                xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
                continue;
            }

            BaseSeqViewModel refViewModel = refEl.refViewModel();
            if(refViewModel == null)
                refViewModel = refEl.assemblyViewModel();

            boolean renderLeftOrientationMarker = i == 0;
            boolean renderRightOrientationMarker = i == mRefViewModel.size() - 1;
            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, mReadViewModel, true, Maps.newHashMap(), refViewModel, renderLeftOrientationMarker, renderRightOrientationMarker);
            renderCount += 1;

            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        if(renderCount <= 1)
            return null;

        String readName = mSupportRead.id();
        ChrBaseRegion alignment = new ChrBaseRegion(mSupportRead.chromosome(), mSupportRead.alignmentStart(), mSupportRead.alignmentEnd());
        ChrBaseRegion mateAlignment = null;
        if(!mSupportRead.mateChromosome().equals(NO_CHROMOSOME_NAME))
            mateAlignment = new ChrBaseRegion(
                    mSupportRead.mateChromosome(), mSupportRead.mateAlignmentStart(), mSupportRead.mateAlignmentEnd());

        String cigarStr = mSupportRead.cigar();
        // TODO(mkcmkc): Include mate cigar string.
        String mateCigarStr = null;
        int insertSize = mSupportRead.insertSize();
        // TODO(mkcmkc): Include orientation string.
        String orientationStr = null;
        int mapQ = mSupportRead.mapQual();
        Integer readNM = null;
        // TODO(mkcmkc): Include this in read info panel.
        String consensusTypeAttribute = null;
        String consensusReadAttribute = null;
        Integer secondMapQ = null;
        Integer secondReadNM = null;

        List<Pair<String, String>> extraInfo = Lists.newArrayList();
        extraInfo.add(Pair.of("Mismatch Info:", mSupportRead.mismatchInfo()));

        String trimCountStr = format("%d,%d", mSupportRead.trimCountStart(), mSupportRead.trimCountEnd());
        extraInfo.add(Pair.of("Trim count:", trimCountStr));

        String mismatchesStr = mSupportRead.cachedRead().mismatches().toString();
        extraInfo.add(Pair.of("Mismatches:", mismatchesStr));

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = rawHtml(svgCanvas.getSVGElement());
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(
                readName, alignment, mateAlignment, cigarStr, mateCigarStr, insertSize, orientationStr, mapQ, readNM, consensusTypeAttribute, consensusReadAttribute, secondMapQ, secondReadNM, extraInfo);
        DomContent containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());

        return containerDiv;
    }
}
