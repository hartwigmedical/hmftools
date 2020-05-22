package com.hartwig.hmftools.linx.fusion.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.fusion.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FIVE_GENE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.THREE_GENE;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.KNOWN;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class RnaFusionAnnotator
{
    private final EnsemblDataCache mGeneTransCache;

    public RnaFusionAnnotator(final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;
    }

    public static RnaExonMatchData findExonMatch(final TranscriptData transData, int rnaPosition)
    {
        RnaExonMatchData exonMatch = new RnaExonMatchData(transData.TransName);

        for (int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exonData = transData.exons().get(i);
            final ExonData nextExonData = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

            if (rnaPosition == exonData.ExonEnd || rnaPosition == exonData.ExonStart)
            {
                // skip matches on the last exon
                if(i == 0 && transData.Strand == -1 && rnaPosition == exonData.ExonStart)
                    return exonMatch;
                else if(i == transData.exons().size() - 1 && transData.Strand == POS_STRAND && rnaPosition == exonData.ExonEnd)
                    return exonMatch;

                // position exactly matches the bounds of an exon
                exonMatch.ExonFound = true;
                exonMatch.BoundaryMatch = true;
                exonMatch.ExonRank = exonData.ExonRank;

                if ((transData.Strand == POS_STRAND) == (rnaPosition == exonData.ExonStart))
                {
                    exonMatch.ExonPhase = exonData.ExonPhase;
                }
                else
                {
                    exonMatch.ExonPhase = exonData.ExonPhaseEnd;
                }
                break;
            }

            if (rnaPosition > exonData.ExonStart && rnaPosition < exonData.ExonEnd)
            {
                // position is within the bounds of an exon
                exonMatch.ExonFound = true;
                exonMatch.ExonRank = exonData.ExonRank;
                exonMatch.ExonPhase = exonData.ExonPhase;
                break;
            }

            if (nextExonData != null && rnaPosition > exonData.ExonEnd && rnaPosition < nextExonData.ExonStart)
            {
                exonMatch.ExonFound = true;

                if (transData.Strand == POS_STRAND)
                {
                    exonMatch.ExonRank = exonData.ExonRank;
                    exonMatch.ExonPhase = exonData.ExonPhaseEnd;
                }
                else
                {
                    exonMatch.ExonRank = nextExonData.ExonRank;
                    exonMatch.ExonPhase = nextExonData.ExonPhaseEnd;
                }

                break;
            }
        }

        return exonMatch;
    }

    public static void checkRnaPhasedTranscripts(RnaFusionData rnaFusion)
    {
        if(rnaFusion.getTransExonData()[FS_UPSTREAM].isEmpty() || rnaFusion.getTransExonData()[FS_DOWNSTREAM].isEmpty())
            return;

        for (RnaExonMatchData exonDataUp : rnaFusion.getTransExonData()[FS_UPSTREAM])
        {
            if(rnaFusion.JunctionTypes[FS_UPSTREAM] == KNOWN && !exonDataUp.BoundaryMatch)
                continue;

            for(RnaExonMatchData exonDataDown : rnaFusion.getTransExonData()[FS_DOWNSTREAM])
            {
                if(rnaFusion.JunctionTypes[FS_DOWNSTREAM] == KNOWN && !exonDataDown.BoundaryMatch)
                    continue;

                boolean phaseMatched = exonDataUp.ExonPhase == exonDataDown.ExonPhase;

                if(phaseMatched && !rnaFusion.hasRnaPhasedFusion())
                {
                    LNX_LOGGER.debug("rnaFusion({}) juncTypes(up={} down={}) transUp({}) transDown({} phase({}) matched({})",
                            rnaFusion.name(), rnaFusion.JunctionTypes[FS_UPSTREAM], rnaFusion.JunctionTypes[FS_DOWNSTREAM],
                            exonDataUp.TransName, exonDataDown.TransName, exonDataUp.ExonPhase, phaseMatched);

                    rnaFusion.setRnaPhasedFusionData(exonDataUp, exonDataDown);
                }
            }
        }
    }

    public static void setReferenceFusionData(final KnownFusionData refFusionData, RnaFusionData rnaFusion)
    {
        if(refFusionData == null)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
            return;
        }

        for(final String[] genePair : refFusionData.knownPairs())
        {
            if (genePair[FIVE_GENE].equals(rnaFusion.GeneNames[FS_UPSTREAM]) && genePair[THREE_GENE].equals(rnaFusion.GeneNames[FS_DOWNSTREAM]))
            {
                rnaFusion.setKnownType(REPORTABLE_TYPE_KNOWN);
                return;
            }
        }

        boolean fivePrimeProm = refFusionData.hasPromiscuousFiveGene(rnaFusion.GeneNames[FS_UPSTREAM]);
        boolean threePrimeProm = refFusionData.hasPromiscuousThreeGene(rnaFusion.GeneNames[FS_DOWNSTREAM]);

        if(fivePrimeProm && threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_BOTH_PROM);
        }
        else if(fivePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_5P_PROM);
        }
        else if(threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_3P_PROM);
        }
        else
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
        }
    }

    public static boolean isViableBreakend(final SvBreakend breakend, int rnaPosition, byte geneStrand, boolean isUpstream)
    {
        boolean requireHigherBreakendPos = isUpstream ? (geneStrand == 1) : (geneStrand == -1);

        int position = breakend.position();

        int offsetMargin = getHomologyMargin(breakend);

        if(requireHigherBreakendPos)
        {
            if(breakend.orientation() != POS_ORIENT)
                return false;

            // factor in any uncertainty around the precise breakend, eg from homology
            return (position + offsetMargin >= rnaPosition);
        }
        else
        {
            if(breakend.orientation() != NEG_ORIENT)
                return false;

            return (position - offsetMargin <= rnaPosition);
        }
    }

    private static int getHomologyMargin(final SvBreakend breakend)
    {
        // the interval offset could be used in place of half the homology but interpretation of the GRIDSS value needs to be understood first
        /*
        int offsetMargin = requireHigherBreakendPos ?
                (breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetEnd() : breakend.getSV().getSvData().endIntervalOffsetEnd())
                : (breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetStart() : breakend.getSV().getSvData().endIntervalOffsetStart());
        */

        int homologyLength = breakend.usesStart() ?
                breakend.getSV().getSvData().startHomologySequence().length() : breakend.getSV().getSvData().endHomologySequence().length();

        return (homologyLength / 2) + (homologyLength % 2);
    }

    public static boolean positionMatch(final SvBreakend breakend, int rnaPosition)
    {
        // checks for a position match within the bounds of uncertainty
        int offsetMargin = getHomologyMargin(breakend);

        return abs(breakend.position() - rnaPosition) <= offsetMargin;
    }

    public boolean isTranscriptBreakendViableForRnaBoundary(
            final Transcript trans, boolean isUpstream, int breakendPosition,
            int rnaPosition, boolean exactRnaPosition)
    {
        // breakend must fall at or before the RNA boundary but not further upstream than the previous splice acceptor

        // if the RNA boundary is at or before the 2nd exon (which has the first splice acceptor), then the breakend can
        // be upstream as far the previous gene or 100K
        final TranscriptData transData = mGeneTransCache.getTranscriptData(trans.gene().StableId, trans.StableId);

        if (transData == null || transData.exons().isEmpty())
            return false;

        int strand = trans.gene().Strand;

        // first find the matching exon boundary for this RNA fusion boundary
        for (int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exonData = transData.exons().get(i);
            final ExonData prevExonData = i > 0 ? transData.exons().get(i - 1) : null;
            final ExonData nextExonData = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

            if (isUpstream)
            {
                // first check if at an exon boundary or before the start of the next exon and after the start of this one
                if(strand == POS_STRAND)
                {
                    if ((rnaPosition == exonData.ExonEnd)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition > exonData.ExonStart && rnaPosition < nextExonData.ExonStart))
                    {
                        // in which case check whether the breakend is before the next exon's splice acceptor
                        if (nextExonData != null)
                        {
                            return breakendPosition < nextExonData.ExonStart;
                        }

                        // can't take the last exon
                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.ExonStart)
                            || (!exactRnaPosition && prevExonData != null && rnaPosition < exonData.ExonEnd && rnaPosition > prevExonData.ExonEnd))
                    {
                        if(prevExonData != null)
                        {
                            return breakendPosition > prevExonData.ExonEnd;
                        }

                        return false;
                    }
                }
            }
            else
            {
                if((strand == POS_STRAND && rnaPosition <= exonData.ExonStart && exonData.ExonRank <= 2)
                || (strand == NEG_STRAND && rnaPosition >= exonData.ExonEnd && exonData.ExonRank <= 2))
                {
                    int breakendDistance = abs(breakendPosition - rnaPosition);

                    if(breakendDistance > PRE_GENE_PROMOTOR_DISTANCE || trans.hasNegativePrevSpliceAcceptorDistance())
                        return false;
                    else
                        return true;
                }

                if(strand == POS_STRAND)
                {
                    if ((rnaPosition == exonData.ExonStart)
                    || (!exactRnaPosition && prevExonData != null && rnaPosition > prevExonData.ExonStart && rnaPosition < exonData.ExonStart))
                    {
                        if(prevExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition > prevExonData.ExonStart;
                        }

                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.ExonEnd)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition < nextExonData.ExonEnd && rnaPosition > exonData.ExonEnd))
                    {
                        if(nextExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition < nextExonData.ExonStart;
                        }

                        return false;
                    }
                }
            }
        }

        return false;
    }


}
