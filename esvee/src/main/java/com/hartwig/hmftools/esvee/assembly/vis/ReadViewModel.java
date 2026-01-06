package com.hartwig.hmftools.esvee.assembly.vis;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.SvgRender.BOX_PADDING;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.READ_HEIGHT_PX;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
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

import org.apache.commons.lang3.tuple.Pair;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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

    public static ReadViewModel create(final Cigar sequenceCigar, final List<Pair<CigarOperator, RefSegmentViewModel>> refViewModel,
            final JunctionAssembly assembly, final SupportRead read)
    {
        Junction junction = assembly.junction();
        SAMRecord record = read.cachedRead().bamRecord();
        int unclippedStart = junction.Position - read.junctionReadStartDistance();
        for(CigarElement cigarEl : record.getCigar().getCigarElements())
        {
            if(cigarEl.getOperator() == I)
                unclippedStart += cigarEl.getLength();
            else if(cigarEl.getOperator() == D)
                unclippedStart -= cigarEl.getLength();
        }

        // find corresponding ref segment and the remote ref segment
        RefSegmentViewModel readSegmentViewModel = null;
        RefSegmentViewModel remoteSegmentViewModel = null;
        int readViewModelIdx = -1;
        int remoteViewModelIdx = -1;
        for(int i = 0; i < refViewModel.size(); i++)
        {
            CigarOperator cigarOp = refViewModel.get(i).getLeft();
            if(cigarOp != M)
                continue;

            RefSegmentViewModel segmentViewModel = refViewModel.get(i).getRight();
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
        List<SegmentViewModel> segmentViewModels = Lists.newArrayList();
        while(segmentViewModels.size() < refViewModel.size())
            segmentViewModels.add(null);

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

        List<BaseViewModel> insertBaseViewModels = null;
        if(remoteBaseViewModels.size() <= insertLength)
        {
            if(!remoteBaseViewModels.isEmpty())
                insertBaseViewModels = remoteBaseViewModels;

            remoteBaseViewModels = Lists.newArrayList();
        }
        else
        {
            if(junction.Orient == FORWARD)
            {
                if(insertLength > 0)
                    insertBaseViewModels = remoteBaseViewModels.subList(0, insertLength);

                remoteBaseViewModels = remoteBaseViewModels.subList(insertLength, remoteBaseViewModels.size());
            }
            else
            {
                if(insertLength > 0)
                    insertBaseViewModels = remoteBaseViewModels.subList(
			    remoteBaseViewModels.size() - insertLength, remoteBaseViewModels.size());

                remoteBaseViewModels = remoteBaseViewModels.subList(0, remoteBaseViewModels.size() - insertLength);
            }
        }

        if(insertBaseViewModels != null)
        {
            RefSegmentViewModel insertSegmentViewModel = refViewModel.get(1).getRight();
            int posStart = readViewModelIdx == 0 ? 1 : 1 + insertLength - insertBaseViewModels.size();
            BaseSeqViewModel insertReadViewModel = new BaseSeqViewModel(insertBaseViewModels, null, posStart, null, null);
            segmentViewModels.set(1, new SegmentViewModel(insertSegmentViewModel.viewRegion(), insertReadViewModel, insertSegmentViewModel.viewModel()));
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
        int i = 0;
        for(SegmentViewModel segmentViewModel : mSegmentViewModels)
        {
            BaseRegion viewRegion = segmentViewModel.viewRegion;
            BaseSeqViewModel readViewModel = segmentViewModel.readViewModel;
            BaseSeqViewModel refViewModel = segmentViewModel.refViewModel;

            svgCanvas.setTransform(initTransform);

            boolean renderLeftOrientationMarker = i == 0;
            boolean renderRightOrientationMarker = i == mSegmentViewModels.size() - 1;
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, readViewModel, true, Maps.newHashMap(), refViewModel, renderLeftOrientationMarker, renderRightOrientationMarker);

            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
            i += 1;
        }

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = rawHtml(svgCanvas.getSVGElement());
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(mRead, null);
        DomContent containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());

        return containerDiv;
    }
}
