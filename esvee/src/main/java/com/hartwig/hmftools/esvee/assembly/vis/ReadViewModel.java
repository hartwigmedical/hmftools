package com.hartwig.hmftools.esvee.assembly.vis;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.SvgRender.BOX_PADDING;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.READ_HEIGHT_PX;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static j2html.TagCreator.div;
import static j2html.TagCreator.rawHtml;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisualiser.RefSegmentViewModel;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.CigarElement;
import j2html.tags.DomContent;

public final class ReadViewModel
{
    private final SupportRead mSupportRead;
    private final List<RefSegmentViewModel> mRefViewModel;
    private final BaseSeqViewModel mReadViewModel;

    private ReadViewModel(final SupportRead supportRead, final List<RefSegmentViewModel> refViewModel, final BaseSeqViewModel readViewModel)
    {
        mSupportRead = supportRead;
        mRefViewModel = refViewModel;
        mReadViewModel = readViewModel;
    }

    public static ReadViewModel create(final List<RefSegmentViewModel> refViewModel, final SupportRead read, final JunctionAssembly junctionAssembly)
    {
        boolean readNegativeStrandFlag = read.cachedRead().bamRecord().getReadNegativeStrandFlag();
        byte[] readBases = read.cachedRead().getBases();
        byte[] readBaseQuals = read.cachedRead().getBaseQuality();
        if(read.fullAssemblyOrientation() == REVERSE)
        {
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
        for(RefSegmentViewModel refViewModel : mRefViewModel)
            totalBoxWidth += refViewModel.viewRegion().baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        int renderCount = 0;
        for(int i = 0; i < mRefViewModel.size(); i += 1)
        {
            RefSegmentViewModel refViewModel = mRefViewModel.get(i);
            BaseRegion viewRegion = refViewModel.viewRegion();

            if(!viewRegion.overlaps(readRegion))
            {
                xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
                continue;
            }

            boolean renderLeftOrientationMarker = i == 0;
            boolean renderRightOrientationMarker = i == mRefViewModel.size() - 1;
            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, mReadViewModel, true, Maps.newHashMap(), refViewModel.viewModel(), renderLeftOrientationMarker, renderRightOrientationMarker);
            renderCount += 1;

            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        if(renderCount <= 1)
            return null;

        Map<String, String> extraInfo = Maps.newHashMap();
        extraInfo.put("Mismatch Info:", mSupportRead.mismatchInfo());

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = rawHtml(svgCanvas.getSVGElement());
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(mSupportRead.cachedRead().bamRecord(), null, extraInfo);
        DomContent containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());

        return containerDiv;
    }
}
