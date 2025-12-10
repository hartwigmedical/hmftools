package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Codons.codonToAminoAcid;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.sage.vis.AminoAcidEvent.START_LOST;
import static com.hartwig.hmftools.sage.vis.AminoAcidEvent.STOP;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.GeneRegionViewModel;
import com.hartwig.hmftools.sage.common.SageVariant;

public final class AminoAcidUtil
{
    private AminoAcidUtil() {}

    public static List<GeneRegionViewModel> getGeneRegions(final TranscriptData transcriptExons)
    {
        boolean posStrand = transcriptExons.posStrand();
        List<GeneRegionViewModel> geneRegions = Lists.newArrayList();
        BaseRegion transcriptCodingRegion = new BaseRegion(transcriptExons.CodingStart, transcriptExons.CodingEnd);
        List<ExonData> exons = transcriptExons.exons();
        ExonData prevExon = null;
        int idxStart = posStrand ? 0 : exons.size() - 1;
        int inc = posStrand ? 1 : -1;
        int aaPos = 1;
        int phase = 0;
        for(int i = idxStart; i >= 0 && i < exons.size(); i += inc)
        {
            ExonData exon = exons.get(i);
            if(prevExon != null)
            {
                if(posStrand)
                    geneRegions.add(new GeneRegionViewModel.IntronicRegionViewModel(new BaseRegion(prevExon.End + 1, exon.Start - 1), 0));
                else
                    geneRegions.add(new GeneRegionViewModel.IntronicRegionViewModel(new BaseRegion(exon.End + 1, prevExon.Start - 1), 0));
            }

            prevExon = exon;

            BaseRegion exonRegion = new BaseRegion(exon.Start, exon.End);
            if(!transcriptCodingRegion.overlaps(exonRegion))
            {
                geneRegions.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(exonRegion, 0));
                continue;
            }

            int codingStart = max(transcriptCodingRegion.start(), exonRegion.start());
            int codingEnd = min(transcriptCodingRegion.end(), exonRegion.end());
            if(posStrand)
            {
                if(exon.Start < codingStart)
                    geneRegions.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(exon.Start, codingStart - 1), 0));
            }
            else
            {
                if(exon.End > codingEnd)
                    geneRegions.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(codingEnd + 1, exon.End), 0));
            }

            BaseRegion codingRegion = new BaseRegion(codingStart, codingEnd);
            int pos = posStrand ? codingStart : codingEnd;
            while(codingRegion.containsPosition(pos))
            {
                int nucLength = 3 - phase;
                int end = posStrand
                        ? min(pos + nucLength - 1, codingRegion.end())
                        : max(pos - nucLength + 1, codingRegion.start());
                int len = end >= pos ? end - pos + 1 : pos - end + 1;
                phase = (len + phase) % 3;
                BaseRegion region = new BaseRegion(min(pos, end), max(pos, end));
                geneRegions.add(new GeneRegionViewModel.AminoAcidViewModel(region, aaPos, '.', '.', 0, false));
                if(phase == 0)
                    aaPos++;

                pos = posStrand ? end + 1 : end - 1;
            }

            if(posStrand)
            {
                if(exon.End > codingEnd)
                    geneRegions.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(codingEnd + 1, exon.End), 0));
            }
            else
            {
                if(exon.Start < codingStart)
                    geneRegions.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(exon.Start, codingStart - 1), 0));
            }

        }

        return geneRegions;
    }

    public static List<GeneRegionViewModel> getRefGeneRegionViewModels(final TranscriptData transcriptExons,
            final TranscriptAminoAcids transcriptAminoAcids, final List<GeneRegionViewModel> geneRegions)
    {
        List<GeneRegionViewModel> viewModels = Lists.newArrayList();
        boolean posStrand = transcriptExons.posStrand();
        String aminoAcids = transcriptAminoAcids.AminoAcids;
        for(GeneRegionViewModel geneRegion : geneRegions)
        {
            if(geneRegion instanceof GeneRegionViewModel.IntronicRegionViewModel)
            {
                viewModels.add(geneRegion);
                continue;
            }

            if(geneRegion instanceof GeneRegionViewModel.NonCodingExonicRegionViewModel)
            {
                viewModels.add(geneRegion);
                continue;
            }

            if(geneRegion instanceof GeneRegionViewModel.AminoAcidViewModel aaRegion)
            {
                int aaPos = aaRegion.aminoAcidPos();
                char refAcid = aminoAcids.charAt(aaPos - 1);
                viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(aaRegion.region(), aaPos, refAcid, '.', 0, aaPos == 1));
                continue;
            }

            throw new IllegalArgumentException(format("Cannot handle type: %s", geneRegion.getClass()));
        }

        if(!posStrand)
            Collections.reverse(viewModels);

        return viewModels;
    }

    public static List<GeneRegionViewModel> getAltGeneRegionViewModels(final TranscriptData transcriptExons,
            final TranscriptAminoAcids transcriptAminoAcids, final List<GeneRegionViewModel> geneRegions,
            final List<AminoAcidEvent> events, final BaseRegion renderRegion, final RefGenomeSource refGenome, final SageVariant variant)
    {
        List<GeneRegionViewModel> viewModels = Lists.newArrayList();
        boolean posStrand = transcriptExons.posStrand();
        String aminoAcids = transcriptAminoAcids.AminoAcids;
        List<AminoAcidEvent> sortedEvents = Lists.newArrayList(events);
        Collections.sort(sortedEvents, AminoAcidEvent::compare);
        Deque<AminoAcidEvent> eventQ = new ArrayDeque<>(sortedEvents);
        boolean stopCoding = false;
        for(int i = 0; i < geneRegions.size(); i++)
        {
            GeneRegionViewModel geneRegion = geneRegions.get(i);
            if(stopCoding)
            {
                viewModels.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(geneRegion.region(), 0));
                continue;
            }

            if(geneRegion instanceof GeneRegionViewModel.IntronicRegionViewModel)
            {
                viewModels.add(geneRegion);
                continue;
            }

            if(geneRegion instanceof GeneRegionViewModel.NonCodingExonicRegionViewModel)
            {
                viewModels.add(geneRegion);
                continue;
            }

            if(geneRegion instanceof GeneRegionViewModel.AminoAcidViewModel aaRegion)
            {
                int aaPos = aaRegion.aminoAcidPos();
                char refAcid = aminoAcids.charAt(aaPos - 1);
                char altAcid = refAcid;

                while(!eventQ.isEmpty() && eventQ.peek().region().end() < aaPos)
                    eventQ.poll();

                boolean inDel = false;
                boolean skipRest = false;
                Stack<AminoAcidEvent> replayEvents = new Stack<>();
                while(!eventQ.isEmpty() && eventQ.peek().region().start() <= aaPos)
                {
                    AminoAcidEvent event = eventQ.poll();
                    if(event instanceof AminoAcidEvent.AminoAcidSubstitution e)
                    {
                        replayEvents.push(event);
                        altAcid = e.alt();

                        // if stop codon is gained, then skip over the rest, and covert to non-coding exon models
                        if(altAcid == STOP)
                        {
                            stopCoding = true;
                            break;
                        }
                    }
                    else if(event instanceof AminoAcidEvent.AminoAcidIns)
                    {
                        // we ignore these here, and add inserts back in when we repack the models
                    }
                    else if(event instanceof AminoAcidEvent.AminoAcidDel e)
                    {
                        replayEvents.push(event);
                        int delLen = e.aminoAcidEndPos() - e.aminoAcidStartPos() + 1;
                        int j;
                        for(j = i + 1; j < geneRegions.size() && j < i + delLen; j++)
                        {
                            if(!(geneRegions.get(j) instanceof GeneRegionViewModel.AminoAcidViewModel))
                                break;
                        }

                        int delEndIdx = j - 1;
                        int delStart = posStrand ? aaRegion.region().start() : geneRegions.get(delEndIdx).region().start();
                        int delEnd = posStrand ? geneRegions.get(delEndIdx).region().end() : aaRegion.region().end();
                        viewModels.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(delStart, delEnd)));
                        i = delEndIdx;
                        inDel = true;
                        break;
                    }
                    else if(event instanceof AminoAcidEvent.AminoAcidExt)
                    {
                        skipRest = true;
                        if(posStrand)
                        {
                            if(aaRegion.region().start() > renderRegion.end())
                                break;

                            final int refPosStart = min(aaRegion.region().start(), variant.position() - 1);
                            final int refPosEnd = min(renderRegion.end() + 100, refGenome.getChromosomeLength(variant.chromosome()));
                            String refBases = refGenome.getBaseString(variant.chromosome(), refPosStart, refPosEnd);
                            int variantIdx = variant.position() - refPosStart;
                            String altBases = refBases.substring(0, variantIdx)
                                    + variant.alt()
                                    + refBases.substring(variantIdx + variant.ref().length());

                            int idx = aaRegion.region().start() - refPosStart;
                            int currentAaPos = aaRegion.aminoAcidPos();
                            Queue<GeneRegionViewModel> geneRegionQ = new ArrayDeque<>(geneRegions.subList(i, geneRegions.size()));
                            while(idx + 2 < altBases.length() && idx + 2 < refBases.length())
                            {
                                BaseRegion region = new BaseRegion(idx + refPosStart, idx + refPosStart + 2);
                                if(region.start() > renderRegion.end())
                                    break;

                                while(!geneRegionQ.isEmpty() && geneRegionQ.peek().region().end() < region.start())
                                    geneRegionQ.poll();

                                if(geneRegionQ.isEmpty())
                                    throw new IllegalStateException("geneRegionQ is empty");

                                char refAA = '.';
                                String codon;
                                if(geneRegionQ.peek() instanceof GeneRegionViewModel.AminoAcidViewModel)
                                {
                                    codon = refBases.substring(idx, idx + 3);
                                    refAA = codonToAminoAcid(codon);
                                }

                                codon = altBases.substring(idx, idx + 3);
                                if(isStopCodon(codon))
                                {
                                    viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(
                                            region, currentAaPos, refAA, '*', 0, false));
                                    stopCoding = true;
                                    break;
                                }

                                char altAA = codonToAminoAcid(codon);
                                viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(
                                        region, currentAaPos, refAA, altAA, 0, currentAaPos == 1));

                                idx += 3;
                                currentAaPos++;
                            }

                            if(stopCoding)
                            {
                                int pos = idx + 3;
                                if(pos <= renderRegion.end())
                                    viewModels.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(pos, renderRegion.end()), 0));
                            }
                        }
                        else
                        {
                            if(aaRegion.region().end() < renderRegion.start())
                                break;

                            final int refPosStart = max(renderRegion.start() - 100, 1);
                            final int refPosEnd = max(aaRegion.region().end(), variant.position() + 1);
                            String refBases = refGenome.getBaseString(variant.chromosome(), refPosStart, refPosEnd);
                            int variantIdx = variant.position() - refPosStart;
                            String altBases = refBases.substring(0, variantIdx)
                                    + variant.alt()
                                    + refBases.substring(variantIdx + variant.ref().length());

                            int refIdx = aaRegion.region().end() - refPosStart;
                            int altIdx = refIdx + variant.alt().length() - variant.ref().length();
                            int currentAaPos = aaRegion.aminoAcidPos();
                            Queue<GeneRegionViewModel> geneRegionQ = new ArrayDeque<>(geneRegions.subList(i, geneRegions.size()));
                            while(refIdx - 2 >= 0 && altIdx - 2 >= 0)
                            {
                                BaseRegion refRegion = new BaseRegion(refIdx + refPosStart - 2, refIdx + refPosStart);
                                if(refRegion.end() < renderRegion.start())
                                    break;

                                while(!geneRegionQ.isEmpty() && geneRegionQ.peek().region().start() > refRegion.end())
                                    geneRegionQ.poll();

                                if(geneRegionQ.isEmpty())
                                    throw new IllegalStateException("geneRegionQ is empty");

                                char refAA = '.';
                                String forwardBases;
                                String codon;
                                if(geneRegionQ.peek() instanceof GeneRegionViewModel.AminoAcidViewModel)
                                {
                                    forwardBases = refBases.substring(refIdx - 2, refIdx + 1);
                                    codon = reverseComplementBases(forwardBases);
                                    refAA = codonToAminoAcid(codon);
                                }

                                forwardBases = altBases.substring(altIdx - 2, altIdx + 1);
                                codon = reverseComplementBases(forwardBases);
                                if(isStopCodon(codon))
                                {
                                    viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(
                                            refRegion, currentAaPos, refAA, '*', 0, false));
                                    stopCoding = true;
                                    break;
                                }

                                char altAA = codonToAminoAcid(codon);
                                viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(
                                        refRegion, currentAaPos, refAA, altAA, 0, currentAaPos == 1));

                                refIdx -= 3;
                                altIdx -= 3;
                                currentAaPos++;
                            }

                            if(stopCoding)
                            {
                                int pos = refIdx - 3;
                                if(pos >= renderRegion.start())
                                    viewModels.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(renderRegion.start(), pos), 0));
                            }
                        }

                        break;
                    }
                    else if(event instanceof AminoAcidEvent.AminoAcidStartLost)
                    {
                        stopCoding = true;
                        altAcid = START_LOST;
                        break;
                    }
                    else
                    {
                        throw new RuntimeException(format("Cannot handle event of type %s", event.getClass()));
                    }
                }

                while(!replayEvents.isEmpty())
                    eventQ.addFirst(replayEvents.pop());

                if(skipRest)
                    break;

                if(inDel)
                    continue;

                if(!stopCoding || altAcid == STOP)
                    viewModels.add(new GeneRegionViewModel.AminoAcidViewModel(aaRegion.region(), aaPos, refAcid, altAcid, 0, aaPos == 1));
                else
                    viewModels.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(aaRegion.region(), 0));

                continue;
            }

            throw new IllegalArgumentException(format("Cannot handle type: %s", geneRegion.getClass()));
        }

        viewModels = realignViewModels(posStrand, variant, viewModels);
        if(!posStrand)
            Collections.reverse(viewModels);

        return viewModels;
    }

    private static List<GeneRegionViewModel> realignViewModels(
            boolean posStrand, final SageVariant variant, final List<GeneRegionViewModel> viewModels)
    {

        List<GeneRegionViewModel> realignedViewModels = Lists.newArrayList();
        GeneRegionViewModel prevViewModel = null;
        for(GeneRegionViewModel viewModel : viewModels)
        {
            if(viewModel instanceof GeneRegionViewModel.DelViewModel)
                continue;

            if(prevViewModel == null)
            {
                realignedViewModels.add(viewModel);
            }
            else if(posStrand)
            {
                int len = viewModel.region().baseLength();
                int start = prevViewModel.region().end() + 1;
                int end = start + len - 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
            }
            else
            {
                int len = viewModel.region().baseLength();
                int end = prevViewModel.region().start() - 1;
                int start = end - len + 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
            }

            prevViewModel = realignedViewModels.get(realignedViewModels.size() - 1);
        }

        int indelSize = variant.alt().length() - variant.ref().length();
        if(indelSize == 0)
            return realignedViewModels;

        if(indelSize < 0)
            return realignViewModelsAroundDel(posStrand, variant, realignedViewModels);

        return realignViewModelsAroundIns(posStrand, variant, realignedViewModels);
    }

    private static List<GeneRegionViewModel> realignViewModelsAroundDel(
            boolean posStrand, final SageVariant variant, final List<GeneRegionViewModel> viewModels)
    {
        int delLen = variant.ref().length() - variant.alt().length();
        BaseRegion delRegion = new BaseRegion(variant.position() + 1, variant.position() + delLen);
        List<GeneRegionViewModel> realignedViewModels = Lists.newArrayList();
        boolean delHandled = false;
        GeneRegionViewModel prevViewModel = null;
        for(GeneRegionViewModel viewModel : viewModels)
        {
            if((posStrand && viewModel.region().end() < delRegion.start()) ||
                    (!posStrand && viewModel.region().start() > delRegion.end()))
            {
                realignedViewModels.add(viewModel);
                prevViewModel = viewModel;
                continue;
            }

            if(delHandled)
            {
                if(prevViewModel == null)
                {
                    realignedViewModels.add(viewModel);
                }
                else if(posStrand)
                {
                    int len = viewModel.region().baseLength();
                    int start = prevViewModel.region().end() + 1;
                    int end = start + len - 1;
                    realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
                }
                else
                {
                    int len = viewModel.region().baseLength();
                    int end = prevViewModel.region().start() - 1;
                    int start = end - len + 1;
                    realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
                }

                prevViewModel = realignedViewModels.get(realignedViewModels.size() - 1);
                continue;
            }

            delHandled = true;
            if(posStrand)
            {
                if(delRegion.start() == viewModel.region().start())
                {
                    int start = viewModel.region().start();
                    int end = start + delLen - 1;
                    realignedViewModels.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(start, end)));

                    int len = viewModel.region().baseLength();
                    start = end + 1;
                    end = start + len - 1;
                    realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));

                    prevViewModel = realignedViewModels.get(realignedViewModels.size() - 1);
                    continue;
                }

                int basesBeforeDel = delRegion.start() - viewModel.region().start();
                int basesAfterDel = 3 - basesBeforeDel;
                int start = viewModel.region().start();
                int end = start + basesBeforeDel - 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));

                start = end + 1;
                end = start + delLen - 1;
                realignedViewModels.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(start, end)));

                start = end + 1;
                end = start + basesAfterDel - 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
            }
            else
            {
                if(delRegion.end() == viewModel.region().end())
                {
                    int end = viewModel.region().end();
                    int start = end - delLen + 1;
                    realignedViewModels.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(start, end)));

                    int len = viewModel.region().baseLength();
                    end = start - 1;
                    start = end - len + 1;
                    realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));

                    prevViewModel = realignedViewModels.get(realignedViewModels.size() - 1);
                    continue;
                }

                int basesBeforeDel = viewModel.region().end() - delRegion.end();
                int basesAfterDel = 3 - basesBeforeDel;
                int end = viewModel.region().end();
                int start = end - basesBeforeDel + 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));

                end = start - 1;
                start = end - delLen + 1;
                realignedViewModels.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(start, end)));

                end = start - 1;
                start = end - basesAfterDel + 1;
                realignedViewModels.add(viewModel.updateRegion(new BaseRegion(start, end)));
            }

            prevViewModel = realignedViewModels.get(realignedViewModels.size() - 1);
        }

        return realignedViewModels;
    }

    private static List<GeneRegionViewModel> realignViewModelsAroundIns(
            boolean posStrand, final SageVariant variant, final List<GeneRegionViewModel> viewModels)
    {
        int insAfterPos = variant.position();
        List<GeneRegionViewModel> realignedViewModels = Lists.newArrayList();
        boolean insHandled = false;
        for(GeneRegionViewModel viewModel : viewModels)
        {
            if((posStrand && viewModel.region().end() < insAfterPos) ||
                    (!posStrand && viewModel.region().start() > insAfterPos + 1))
            {
                realignedViewModels.add(viewModel);
                continue;
            }

            if(insHandled)
            {
                realignedViewModels.add(viewModel);
                continue;
            }

            insHandled = true;
            int basesBeforeRightInsert;
            if(posStrand)
                basesBeforeRightInsert = insAfterPos - viewModel.region().start() + 1;
            else
                basesBeforeRightInsert = viewModel.region().end() - insAfterPos;

            realignedViewModels.add(viewModel.updateBasesBeforeRightInsert(basesBeforeRightInsert));
        }

        return realignedViewModels;
    }

    public static String getGeneRegionLabel(final TranscriptData transcriptExons, int position)
    {
        boolean posStrand = transcriptExons.posStrand();
        List<ExonData> exons = transcriptExons.exons();
        int exonStart = posStrand ? 0 : exons.size() - 1;
        int exonEnd = posStrand ? exons.size() - 1 : 0;
        int exonInc = posStrand ? 1 : -1;
        int exonIdx = 0;
        boolean inExon = false;
        for(int i = exonStart; i >= 0 && i < exons.size(); i += exonInc)
        {
            ExonData exon = exons.get(i);
            if(posStrand && position < exon.Start)
                break;

            if(!posStrand && exon.End < position)
                break;

            exonIdx++;
            if(posStrand && position <= exon.End)
            {
                inExon = true;
                break;
            }

            if(!posStrand && exon.Start <= position)
            {
                inExon = true;
                break;
            }

            if(i == exonEnd)
                return "Downstream";
        }

        if(exonIdx == 0)
            return "Upstream";

        if(inExon)
            return "exon " + exonIdx;

        return "intron " + exonIdx;
    }
}
