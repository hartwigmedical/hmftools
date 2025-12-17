package com.hartwig.hmftools.esvee.assembly.vis;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.SvgRender.BOX_PADDING;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.READ_HEIGHT_PX;

import static htsjdk.samtools.CigarOperator.I;
import static j2html.TagCreator.div;
import static j2html.TagCreator.rawHtml;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.BaseViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisualiser.RefSegmentViewModel;

import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import j2html.tags.DomContent;

public final class ReadViewModel
{
    private record SegmentViewModel(BaseRegion viewRegion, BaseSeqViewModel readViewModel, BaseSeqViewModel refViewModel) {}

    private final SAMRecord mRead;
    private final List<SegmentViewModel> mSegmentViewModels;

    private ReadViewModel(final SAMRecord read, final List<SegmentViewModel> segmentViewModels)
    {
        mRead = read;
        mSegmentViewModels = segmentViewModels;
    }

    public static ReadViewModel create(final Cigar sequenceCigar, final List<RefSegmentViewModel> refViewModel,
            final JunctionAssembly assembly, final SupportRead read)
    {
        Junction junction = assembly.junction();
        int unclippedStart = junction.Position - read.junctionReadStartDistance();
        SAMRecord record = read.cachedRead().bamRecord();

        // find corresponding ref segment and the remote ref segment
        RefSegmentViewModel readSegmentViewModel = null;
        RefSegmentViewModel remoteSegmentViewModel = null;
        int readViewModelIdx = -1;
        int remoteViewModelIdx = -1;
        for(int i = 0; i < refViewModel.size(); i++)
        {
            RefSegmentViewModel segmentViewModel = refViewModel.get(i);
            Junction segmentJunction = segmentViewModel.assembly().junction();

            boolean matches = junction.Chromosome.equals(segmentJunction.Chromosome);
            if(junction.Position != segmentJunction.Position)
                matches = false;

            if(junction.Orient != segmentJunction.Orient)
                matches = false;

            if(matches)
            {
                readViewModelIdx = i;
                readSegmentViewModel = segmentViewModel;
                continue;
            }

            remoteViewModelIdx = i;
            remoteSegmentViewModel = segmentViewModel;
        }

        // construct segment view models
        List<SegmentViewModel> segmentViewModels = Lists.newArrayList(null, null);

        BaseSeqViewModel readViewModel = BaseSeqViewModel.fromRead(record, unclippedStart);
        segmentViewModels.set(readViewModelIdx, new SegmentViewModel(readSegmentViewModel.viewRegion(), readViewModel, readSegmentViewModel.viewModel()));

        int remoteStart = junction.Orient == FORWARD ? junction.Position + 1 : readViewModel.FirstBasePos;
        int remoteEnd = junction.Orient == FORWARD ? readViewModel.LastBasePos : junction.Position - 1;
        List<BaseViewModel> remoteBaseViewModels = Lists.newArrayList();
        if(junction.Orient == FORWARD)
        {
            BaseViewModel baseViewModel = readViewModel.getBase(junction.Position);
            if(baseViewModel.hasCharBase())
            {
                List<Character> rightInsertBases = baseViewModel.rightInsertBases();
                List<Integer> rightInsertBaseQs = baseViewModel.rightInsertBaseQs();
                for(int i = 0; i < rightInsertBases.size(); i++)
                    remoteBaseViewModels.add(new BaseViewModel(rightInsertBases.get(i), rightInsertBaseQs.get(i)));

                baseViewModel.resetRightInsert();
            }
        }

        for(int pos = remoteStart; pos <= remoteEnd; pos++)
        {
            BaseViewModel baseViewModel = readViewModel.getBase(pos);
            if(!baseViewModel.hasCharBase())
                continue;

            remoteBaseViewModels.add(new BaseViewModel(baseViewModel.charBase(), baseViewModel.baseQ()));
            List<Character> rightInsertBases = baseViewModel.rightInsertBases();
            List<Integer> rightInsertBaseQs = baseViewModel.rightInsertBaseQs();
            for(int i = 0; i < rightInsertBases.size(); i++)
                remoteBaseViewModels.add(new BaseViewModel(rightInsertBases.get(i), rightInsertBaseQs.get(i)));
        }

        int insertLength = 0;
        if(sequenceCigar.getCigarElements().size() == 3)
        {
            if(sequenceCigar.getCigarElements().get(1).getOperator() == I)
                insertLength = sequenceCigar.getCigarElements().get(1).getLength();
        }

        if(remoteBaseViewModels.size() <= insertLength)
        {
            remoteBaseViewModels.clear();
        }
        else
        {
            if(junction.Orient == FORWARD)
                remoteBaseViewModels = remoteBaseViewModels.subList(insertLength, remoteBaseViewModels.size());
            else
                remoteBaseViewModels = remoteBaseViewModels.subList(0, remoteBaseViewModels.size() - insertLength);
        }

        Junction remoteJunction = remoteSegmentViewModel.assembly().junction();
        BaseSeqViewModel remoteReadViewModel;
        if(remoteJunction.Orient == REVERSE)
        {
            remoteReadViewModel = new BaseSeqViewModel(
                    remoteBaseViewModels, remoteJunction.Chromosome, remoteJunction.Position, readViewModel.LeftIsForwardStrand, readViewModel.RightIsForwardStrand);
        }
        else
        {
            int posStart = remoteJunction.Position - remoteBaseViewModels.size() + 1;
            remoteReadViewModel = new BaseSeqViewModel(
                    remoteBaseViewModels, remoteJunction.Chromosome, posStart, readViewModel.LeftIsForwardStrand, readViewModel.RightIsForwardStrand);
        }

        segmentViewModels.set(remoteViewModelIdx, new SegmentViewModel(remoteSegmentViewModel.viewRegion(), remoteReadViewModel, remoteSegmentViewModel.viewModel()));
        return new ReadViewModel(record, segmentViewModels);
    }

    public DomContent render()
    {
        int totalBoxWidth = 0;
        for(SegmentViewModel segmentViewModel : mSegmentViewModels)
            totalBoxWidth += segmentViewModel.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(SegmentViewModel segmentViewModel : mSegmentViewModels)
        {
            BaseRegion viewRegion = segmentViewModel.viewRegion;
            BaseSeqViewModel readViewModel = segmentViewModel.readViewModel;
            BaseSeqViewModel refViewModel = segmentViewModel.refViewModel;

            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, readViewModel, true, Maps.newHashMap(), refViewModel);

            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = rawHtml(svgCanvas.getSVGElement());
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(mRead, null);
        DomContent containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());

        return containerDiv;
    }
}
