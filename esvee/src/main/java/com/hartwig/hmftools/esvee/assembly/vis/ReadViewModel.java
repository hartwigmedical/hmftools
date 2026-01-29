package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.collapseCigarOps;
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
import java.util.Collections;
import java.util.List;
import java.util.NavigableMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffInfo;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffType;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisualiser.PairedSegmentViewModel;
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
    private final PairedSegmentViewModel mRefViewModel;
    private final BaseSeqViewModel mReadViewModel;
    private final int mIndelOffset;
    private final boolean mIsReverseComplemented;
    private final int mRefViewModelIndex;
    private final String mFullAssemblyStr;

    private ReadViewModel(final SupportRead supportRead, final PairedSegmentViewModel refViewModel, final BaseSeqViewModel readViewModel, int indelOffset, boolean isReverseComplemented, int refViewModelIndex, final String fullAssemblyStr)
    {
        mSupportRead = supportRead;
        mRefViewModel = refViewModel;
        mReadViewModel = readViewModel;
        mIndelOffset = indelOffset;
        mIsReverseComplemented = isReverseComplemented;
        mRefViewModelIndex = refViewModelIndex;
        mFullAssemblyStr = fullAssemblyStr;
    }

    public int firstBaseIndex()
    {
        return mReadViewModel.FirstBasePos;
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

    private static Integer getRefViewModelIndex(final List<SegmentViewModel> refViewModel, final SupportRead read)
    {
        BaseRegion refAlignment = new BaseRegion(read.alignmentStart(), read.alignmentEnd());
        BaseRegion assemblyAlignment = new BaseRegion(read.fullAssemblyIndexStart(), read.fullAssemblyIndexEnd());
        for(int i = 0; i < refViewModel.size(); i++)
        {
            SegmentViewModel refEl = refViewModel.get(i);
            String chromosome = refEl.chromosome();
            BaseRegion refRegion = refEl.refRegion();
            BaseRegion assemblyRegion = refEl.assemblyRegion();
            if(refRegion == null)
                continue;

            if(read.chromosome().equals(chromosome) && refAlignment.overlaps(refRegion) && assemblyAlignment.overlaps(assemblyRegion))
                return i;
        }

        return null;
    }

    private record CorrectIndelsResult(List<CigarElement> cigarEls, int indexOffset) {}

    private static CorrectIndelsResult correctIndels(final List<SegmentViewModel> refViewModel, final SupportRead read, final List<CigarElement> cigarEls)
    {

        Integer refViewModelIndex = getRefViewModelIndex(refViewModel, read);
        if(refViewModelIndex == null)
        {
            SV_LOGGER.error("Cannot find matching junction for read: {}", read.id());
            System.exit(0);
        }

        boolean isBuiltForward = refViewModelIndex == 0;
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
                        int baseIdx = mismatch.ReadIndex + (isBuiltForward ? 0 : 1) + (i + 1) * (isBuiltForward ? 1 : -1);
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

    private record HandleDelResult(BaseSeqViewModel readViewModel, int indelOffset, List<CigarElement> cigarEls) {}

    @Nullable
    private static HandleDelResult handleDel(final PairedSegmentViewModel refViewModel, final SupportRead read, int indelOffset, final List<CigarElement> cigarEls, final byte[] readBases, final byte[] readBaseQuals, boolean readNegativeStrandFlag)
    {
        int delLength = refViewModel.viewModels().get(0).leftDelLength() + refViewModel.viewModels().get(refViewModel.viewModels().size() - 1).leftDelLength();
        if(delLength <= 0)
            return null;

        List<CigarOperator> cigarOps = Lists.newArrayList();
        for(CigarElement cigarEl : cigarEls)
        {
            for(int i = 0; i < cigarEl.getLength(); i++)
                cigarOps.add(cigarEl.getOperator());
        }

        int readStartIdx = read.fullAssemblyIndexStart() + indelOffset + refViewModel.readStartOffset();
        int lastIdxBeforeBreak = refViewModel.viewModels().get(0).assemblyViewModel().LastBasePos;
        int opIdx = 0;
        List<CigarOperator> newCigarOps = Lists.newArrayList();
        int baseIdx = readStartIdx;
        if(baseIdx <= lastIdxBeforeBreak + delLength)
        {
            while(baseIdx <= lastIdxBeforeBreak && opIdx < cigarOps.size())
            {
                CigarOperator currentCigarOp = cigarOps.get(opIdx);
                opIdx++;
                newCigarOps.add(currentCigarOp);
                if(currentCigarOp.consumesReadBases())
                    baseIdx++;
            }
        }

        if(baseIdx <= lastIdxBeforeBreak)
            return null;

        int readDelLength = delLength - (baseIdx - lastIdxBeforeBreak - 1);
        for(int i = 0; i < readDelLength; i++)
            newCigarOps.add(D);

        while(opIdx < cigarOps.size())
        {
            newCigarOps.add(cigarOps.get(opIdx));
            opIdx++;
        }

        List<CigarElement> newCigarEls = collapseCigarOps(newCigarOps);

        BaseSeqViewModel unshiftedReadViewModel = BaseSeqViewModel.create(
                read.fullAssemblyIndexStart() + indelOffset + refViewModel.readStartOffset(), newCigarEls, readBases, readBaseQuals, readNegativeStrandFlag);
        unshiftedReadViewModel.clearSoftClips();
        int unshiftedMismatchCount = 0;
        for(SegmentViewModel refEl : refViewModel.viewModels())
        {
            if(refEl.refViewModel() != null)
                unshiftedMismatchCount += unshiftedReadViewModel.mismatchCount(refEl.refViewModel());
        }

        indelOffset -= readDelLength;
        BaseSeqViewModel shiftedReadViewModel = BaseSeqViewModel.create(
                read.fullAssemblyIndexStart() + indelOffset + refViewModel.readStartOffset(), newCigarEls, readBases, readBaseQuals, readNegativeStrandFlag);
        shiftedReadViewModel.clearSoftClips();
        int shiftedMismatchCount = 0;
        for(SegmentViewModel refEl : refViewModel.viewModels())
        {
            if(refEl.refViewModel() != null)
                shiftedMismatchCount += shiftedReadViewModel.mismatchCount(refEl.refViewModel());
        }

        if(unshiftedMismatchCount <= shiftedMismatchCount)
        {
            indelOffset += readDelLength;
            return new HandleDelResult(unshiftedReadViewModel, indelOffset, newCigarEls);
        }

        return new HandleDelResult(shiftedReadViewModel, indelOffset, newCigarEls);
    }

    @Nullable
    public static ReadViewModel create(final PairedSegmentViewModel refViewModel, final SupportRead read, final JunctionAssembly junctionAssembly, final String fullAssemblyStr)
    {
        boolean readNegativeStrandFlag = read.orientation() == REVERSE;
        byte[] readBases = read.cachedRead().getBases();
        byte[] readBaseQuals = read.cachedRead().getBaseQuality();
        List<CigarElement> cigarEls = read.cachedRead().cigarElements();

        Integer refViewModelIndex = getRefViewModelIndex(refViewModel.viewModels(), read);
        if(refViewModelIndex == null)
            return null;

        SegmentViewModel readRefViewModel = refViewModel.viewModels().get(refViewModelIndex);
        boolean reverseComplemented = readRefViewModel.isRefReversed();
        if(readRefViewModel.isAssemblyReversed())
            reverseComplemented = !reverseComplemented;

        boolean isForwardJunction = junctionAssembly.isForwardJunction();
        if(reverseComplemented)
            isForwardJunction = !isForwardJunction;

        int indelOffset = 0;
        if(isForwardJunction)
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
            CorrectIndelsResult correctIndelsResult = correctIndels(refViewModel.viewModels(), read, cigarEls);
            cigarEls = correctIndelsResult.cigarEls;
            indelOffset += correctIndelsResult.indexOffset;
        }

        if(reverseComplemented)
        {
            readNegativeStrandFlag = !readNegativeStrandFlag;
            readBases = reverseComplementBases(readBases);
            Collections.reverse(cigarEls);

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

        HandleDelResult handleDelResult = handleDel(
                refViewModel, read, indelOffset, cigarEls, readBases, readBases, readNegativeStrandFlag);
        BaseSeqViewModel readViewModel;
        if(handleDelResult == null)
        {
            readViewModel = BaseSeqViewModel.create(
                    read.fullAssemblyIndexStart() + indelOffset + refViewModel.readStartOffset(), cigarEls, readBases, readBaseQuals, readNegativeStrandFlag);
            readViewModel.clearSoftClips();
        }
        else
        {
            readViewModel = handleDelResult.readViewModel;
            indelOffset = handleDelResult.indelOffset;
        }

        return new ReadViewModel(read, refViewModel, readViewModel, indelOffset, reverseComplemented, refViewModelIndex, fullAssemblyStr);
    }

    @Nullable
    public DomContent render()
    {
        BaseRegion readRegion = new BaseRegion(mReadViewModel.FirstBasePos, mReadViewModel.LastBasePos);

        int totalBoxWidth = mRefViewModel.prevBoxWidth();
        for(SegmentViewModel refViewModel : mRefViewModel.viewModels())
            totalBoxWidth += refViewModel.viewRegion().baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = mRefViewModel.prevBoxWidth();
        int renderCount = 0;
        for(int i = 0; i < mRefViewModel.viewModels().size(); i += 1)
        {
            SegmentViewModel refEl = mRefViewModel.viewModels().get(i);
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
            boolean renderRightOrientationMarker = i == mRefViewModel.viewModels().size() - 1;
            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, mReadViewModel, true, Maps.newHashMap(), refViewModel, renderLeftOrientationMarker, renderRightOrientationMarker, refEl.isRefReversed());
            renderCount += 1;

            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        if(renderCount <= 0)
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

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Full assembly index start:", String.valueOf(mSupportRead.fullAssemblyIndexStart())));

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Full assembly seq start:", mFullAssemblyStr.substring(mSupportRead.fullAssemblyIndexStart(), min(mSupportRead.fullAssemblyIndexEnd() + 1, mSupportRead.fullAssemblyIndexStart() + 10))));

        // TODO(mkcmkc): Remove
        SupplementaryReadData suppData = mSupportRead.supplementaryData();
        extraInfo.add(Pair.of("Supplementary Data:", suppData == null ? "None" : suppData.toString()));

        // TODO(mkcmkc): Remove
        String mismatchesStr = mSupportRead.cachedRead().mismatches().toString();
        extraInfo.add(Pair.of("Mismatches:", mismatchesStr));

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Full assembly orientation:", mSupportRead.fullAssemblyOrientation().name()));

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Indel offset:", String.valueOf(mIndelOffset)));

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Reverse complemented:", String.valueOf(mIsReverseComplemented)));

        // TODO(mkcmkc): Remove
        extraInfo.add(Pair.of("Ref view model index:", String.valueOf(mRefViewModelIndex)));

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = rawHtml(svgCanvas.getSVGElement());
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(
                readName, alignment, mateAlignment, cigarStr, mateCigarStr, insertSize, orientationStr, mapQ, readNM, consensusTypeAttribute, consensusReadAttribute, secondMapQ, secondReadNM, extraInfo);
        DomContent containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());

        return containerDiv;
    }
}
