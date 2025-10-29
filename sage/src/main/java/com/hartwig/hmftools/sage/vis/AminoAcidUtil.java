package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.cider.IgTcrGeneFile.Column.posEnd;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent.AminoAcidDel;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent.AminoAcidExt;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent.AminoAcidIns;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent.AminoAcidStartLost;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent.AminoAcidSubstitution;
import com.hartwig.hmftools.sage.vis.GeneRegionViewModel.AminoAcidViewModel;
import com.hartwig.hmftools.sage.vis.GeneRegionViewModel.IntronicRegionViewModel;
import com.hartwig.hmftools.sage.vis.GeneRegionViewModel.NonCodingExonicRegionViewModel;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public final class AminoAcidUtil
{
    private AminoAcidUtil() {}

    private enum GeneRegionType
    {
        NON_CODING_EXON,
        CODING_EXON,
        INTRON
    }

    public static List<Pair<GeneRegionType, BaseRegion>> getGeneRegions(final TranscriptData transcriptExons)
    {
        List<Pair<GeneRegionType, BaseRegion>> exonRegions = Lists.newArrayList();
        BaseRegion transcriptCodingRegion = new BaseRegion(transcriptExons.CodingStart, transcriptExons.CodingEnd);
        List<ExonData> exons = transcriptExons.exons();
        ExonData prevExon = null;
        for(ExonData exon : exons)
        {
            if(prevExon != null)
                exonRegions.add(Pair.of(GeneRegionType.INTRON, new BaseRegion(prevExon.End + 1, exon.Start - 1)));

            prevExon = exon;

            BaseRegion exonRegion = new BaseRegion(exon.Start, exon.End);
            if(!transcriptCodingRegion.overlaps(exonRegion))
            {
                exonRegions.add(Pair.of(GeneRegionType.NON_CODING_EXON, exonRegion));
                continue;
            }

            int codingStart = max(transcriptCodingRegion.start(), exonRegion.start());
            int codingEnd = min(transcriptCodingRegion.end(), exonRegion.end());
            if(exon.Start < codingStart)
                exonRegions.add(Pair.of(GeneRegionType.NON_CODING_EXON, new BaseRegion(exon.Start, codingStart - 1)));

            BaseRegion codingRegion = new BaseRegion(codingStart, codingEnd);
            exonRegions.add(Pair.of(GeneRegionType.CODING_EXON, codingRegion));

            if(exon.End > codingEnd)
                exonRegions.add(Pair.of(GeneRegionType.NON_CODING_EXON, new BaseRegion(codingEnd + 1, exon.End)));
        }

        return exonRegions;
    }

    public static List<GeneRegionViewModel> getGeneRegionViewModels(
            final TranscriptData transcriptExons, final TranscriptAminoAcids transcriptAminoAcids, final List<AminoAcidEvent> events_,
            final BaseSeqViewModel refNucs)
    {
        List<GeneRegionViewModel> viewModels = Lists.newArrayList();
        boolean posStrand = transcriptExons.Strand == (byte) 1;
        String aminoAcids = transcriptAminoAcids.AminoAcids;
        List<Pair<GeneRegionType, BaseRegion>> geneRegions = getGeneRegions(transcriptExons);
        List<AminoAcidEvent> sortedEvents = Lists.newArrayList(events_);
        Collections.sort(sortedEvents, AminoAcidEvent::compare);
        if(!posStrand)
        {
            Collections.reverse(geneRegions);
            Collections.reverse(sortedEvents);
        }

        Queue<AminoAcidEvent> eventQ = new ArrayDeque<>(sortedEvents);
        int nucPhase = 0;
        int aaIdx = 0;
        boolean stopCoding = false;
        for(int i = 0; i < geneRegions.size(); i++)
        {
            Pair<GeneRegionType, BaseRegion> geneRegion = geneRegions.get(i);
            if(geneRegion.getLeft() == GeneRegionType.INTRON)
            {
                viewModels.add(new IntronicRegionViewModel(geneRegion.getRight()));
                continue;
            }

            if(stopCoding || geneRegion.getLeft() == GeneRegionType.NON_CODING_EXON)
            {
                viewModels.add(new NonCodingExonicRegionViewModel(geneRegion.getRight()));
                continue;
            }

            int pos = posStrand ? geneRegion.getRight().start() : geneRegion.getRight().end();
            while(geneRegion.getRight().containsPosition(pos))
            {
                int nucLength = 3 - nucPhase;
                int end = posStrand
                        ? min(pos + nucLength - 1, geneRegion.getRight().end())
                        : max(pos - nucLength + 1, geneRegion.getRight().start());
                int len = end >= pos ? end - pos + 1 : pos - end + 1;
                nucPhase = (len + nucPhase) % 3;
                char refAcid = aminoAcids.charAt(aaIdx);
                int aaPos = aaIdx + 1;
                BaseRegion aaRegion = new BaseRegion(min(pos, end), max(pos, end));
                char altAcid = refAcid;

                while(!eventQ.isEmpty() && eventQ.peek().region().end() < aaPos)
                    eventQ.poll();

                int insertLength = 0;
                while(!eventQ.isEmpty() && eventQ.peek().region().start() <= aaPos)
                {
                    AminoAcidEvent event = eventQ.poll();
                    if(event instanceof AminoAcidSubstitution)
                    {
                        AminoAcidSubstitution e = (AminoAcidSubstitution) event;
                        altAcid = e.alt();

                        // if stop codon is gained, then skip over the rest, and covert to non-coding exon models
                        if(altAcid == '*')
                        {
                            stopCoding = true;
                            break;
                        }
                    }
                    else if(event instanceof AminoAcidIns)
                    {
                        AminoAcidIns e = (AminoAcidIns) event;
                        insertLength += e.length();
                    }
                    else if(event instanceof AminoAcidDel)
                    {
                        altAcid = '.';
                    }
                    else if(event instanceof AminoAcidExt)
                    {
                        // TODO:
                        throw new NotImplementedException("TODO");
                    }
                    else if(event instanceof AminoAcidStartLost)
                    {
                        stopCoding = true;
                        altAcid = '?';
                        break;
                    }
                    else
                    {
                        throw new RuntimeException(format("Cannot handle event of type %s", event.getClass()));
                    }
                }

                if(!stopCoding || altAcid == '*')
                    viewModels.add(new AminoAcidViewModel(aaRegion, aaPos, refAcid, altAcid, insertLength));

                if(stopCoding)
                {
                    int posStart = altAcid == '*' ? (posStrand ? end + 1 : end - 1) : pos;
                    int posEnd = posStrand ? geneRegion.getRight().end() : geneRegion.getRight().start();

                    viewModels.add(new NonCodingExonicRegionViewModel(new BaseRegion(min(posStart, posEnd), max(posStart, posEnd))));
                    break;
                }

                if(nucPhase == 0)
                    aaIdx++;

                pos = posStrand ? end + 1 : end - 1;
            }
        }

        if(!posStrand)
            Collections.reverse(viewModels);

        return viewModels;
    }
}
