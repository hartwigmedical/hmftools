package com.hartwig.hmftools.linx.fusion.rna;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_BOTH;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.KNOWN;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
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

            if (rnaPosition == exonData.End || rnaPosition == exonData.Start)
            {
                // skip matches on the last exon
                if(i == 0 && transData.Strand == -1 && rnaPosition == exonData.Start)
                    return exonMatch;
                else if(i == transData.exons().size() - 1 && transData.Strand == POS_STRAND && rnaPosition == exonData.End)
                    return exonMatch;

                // position exactly matches the bounds of an exon
                exonMatch.ExonFound = true;
                exonMatch.BoundaryMatch = true;
                exonMatch.ExonRank = exonData.Rank;

                if ((transData.Strand == POS_STRAND) == (rnaPosition == exonData.Start))
                {
                    exonMatch.ExonPhase = exonData.PhaseStart;
                }
                else
                {
                    exonMatch.ExonPhase = exonData.PhaseEnd;
                }
                break;
            }

            if (rnaPosition > exonData.Start && rnaPosition < exonData.End)
            {
                // position is within the bounds of an exon
                exonMatch.ExonFound = true;
                exonMatch.ExonRank = exonData.Rank;
                exonMatch.ExonPhase = exonData.PhaseStart;
                break;
            }

            if (nextExonData != null && rnaPosition > exonData.End && rnaPosition < nextExonData.Start)
            {
                exonMatch.ExonFound = true;

                if (transData.Strand == POS_STRAND)
                {
                    exonMatch.ExonRank = exonData.Rank;
                    exonMatch.ExonPhase = exonData.PhaseEnd;
                }
                else
                {
                    exonMatch.ExonRank = nextExonData.Rank;
                    exonMatch.ExonPhase = nextExonData.PhaseEnd;
                }

                break;
            }
        }

        return exonMatch;
    }

    public static void checkRnaPhasedTranscripts(RnaFusionData rnaFusion)
    {
        if(rnaFusion.getTransExonData(FS_UP).isEmpty() || rnaFusion.getTransExonData(FS_DOWN).isEmpty())
            return;

        for (RnaExonMatchData exonDataUp : rnaFusion.getTransExonData(FS_UP))
        {
            if(rnaFusion.JunctionTypes[FS_UP] == KNOWN && !exonDataUp.BoundaryMatch)
                continue;

            for(RnaExonMatchData exonDataDown : rnaFusion.getTransExonData(FS_DOWN))
            {
                if(rnaFusion.JunctionTypes[FS_DOWN] == KNOWN && !exonDataDown.BoundaryMatch)
                    continue;

                boolean phaseMatched = exonDataUp.ExonPhase == exonDataDown.ExonPhase;

                if(phaseMatched && !rnaFusion.hasRnaPhasedFusion())
                {
                    LNX_LOGGER.debug("rnaFusion({}) juncTypes(up={} down={}) transUp({}) transDown({} phase({}) matched({})",
                            rnaFusion.name(), rnaFusion.JunctionTypes[FS_UP], rnaFusion.JunctionTypes[FS_DOWN],
                            exonDataUp.TransName, exonDataDown.TransName, exonDataUp.ExonPhase, phaseMatched);

                    rnaFusion.setRnaPhasedFusionData(exonDataUp, exonDataDown);
                }
            }
        }
    }

    public void correctGeneNames(final KnownFusionCache refFusionData, RnaFusionData rnaFusion)
    {
        // take the first gene where multiple are listed unless one features in the known or promiscuous lists
        final String[] geneNames = rnaFusion.GeneNames;

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(!rnaFusion.GeneIds[fs].isEmpty())
                continue;

            if(rnaFusion.GeneNames[fs].isEmpty())
                continue;

            if(geneNames[fs].contains(";"))
            {
                //	Arriba example: RP11-123H22.1(13969);ATP6V1G1P7(2966)
                String[] genes = geneNames[fs].split(";");
                String selectGene = "";

                for(int i = 0; i < genes.length; ++i)
                {
                    String geneName = genes[i];
                    geneName = geneName.replaceAll("\\([0-9]*\\)", "");

                    if(anyReferenceGeneMatch(refFusionData, geneName))
                    {
                        selectGene = geneName;
                    }
                    else if(selectGene.isEmpty())
                    {
                        selectGene = geneName;
                    }
                }

                geneNames[fs] = selectGene;
            }

            // check that gene names match Ensembl - applicable for StarFusion
            geneNames[fs] = checkAlternateGeneName(geneNames[fs]);

            final GeneData geneData = mGeneTransCache.getGeneDataByName(geneNames[fs]);

            if(geneData == null)
            {
                LNX_LOGGER.warn("geneName({}) gene not found", geneNames[fs]);
                geneNames[fs] = "";
            }
            else
            {
                rnaFusion.GeneIds[fs] = geneData.GeneId;
            }
        }
    }

    private boolean anyReferenceGeneMatch(final KnownFusionCache refFusionData, final String geneName)
    {
        if(refFusionData.getDataByType(KNOWN_PAIR).stream().anyMatch(x -> x.FiveGene.equals(geneName) || x.ThreeGene.equals(geneName)))
            return true;

        return refFusionData.hasPromiscuousFiveGene(geneName) || refFusionData.hasPromiscuousThreeGene(geneName);
    }

    private static String checkAlternateGeneName(final String geneName)
    {
        if(geneName.equals("AC005152.2"))
            return "SOX9-AS1";

        if(geneName.equals("AC016683.6"))
            return "PAX8-AS1";

        if(geneName.equals("AC007092.1"))
            return "LINC01122";

        if(geneName.toUpperCase().equals("C10ORF112"))
            return "MALRD1";

        if(geneName.equals("C5orf50"))
            return "SMIM23";

        if(geneName.equals("CTC-349C3.1"))
            return "C5orf66";

        if(geneName.equals("RP11-500G22.2"))
            return "ATE1-AS1";

        if(geneName.equals("AC018865.8"))
            return "FAR2P1";

        if(geneName.equals("C10orf68"))
            return geneName.toUpperCase();

        if(geneName.equals("C17orf76-AS1"))
            return "FAM211A-AS1";

        if(geneName.equals("IGH@") || geneName.equals("IGH-@"))
            return "IGHJ6";

        if(geneName.equals("IGL@") || geneName.equals("IGL-@"))
            return "IGLC6";

        if(geneName.equals("MKLN1-AS1"))
            return "LINC-PINT";

        if(geneName.equals("PHF15"))
            return "JADE2";

        if(geneName.equals("PHF17"))
            return "JADE1";

        if(geneName.equals("RP11-134P9.1"))
            return "LINC01136";

        if(geneName.equals("RP11-973F15.1"))
            return "LINC01151";

        if(geneName.equals("RP11-115K3.2"))
            return "YWHAEP7";

        if(geneName.equals("RP11-3B12.1"))
            return "POT1-AS1";

        if(geneName.equals("RP11-199O14.1"))
            return "CASC20";

        if(geneName.equals("RP11-264F23.3"))
            return "CCND2-AS1";

        if(geneName.equals("RP11-93L9.1"))
            return "LINC01091";

        if(geneName.equals("AC129929.5"))
            return "CD81-AS1";

        if(geneName.equals("GOLGA2B"))
            return "GOLGA2P5";

        return geneName;
    }

    public static void setReferenceFusionData(final KnownFusionCache refFusionData, RnaFusionData rnaFusion)
    {
        if(refFusionData == null)
        {
            rnaFusion.setKnownType(KnownFusionType.NONE.toString());
            return;
        }

        if(refFusionData.hasKnownFusion(rnaFusion.GeneNames[FS_UP], rnaFusion.GeneNames[FS_DOWN]))
        {
            rnaFusion.setKnownType(KNOWN_PAIR.toString());
            return;
        }

        boolean fivePrimeProm = refFusionData.hasPromiscuousFiveGene(rnaFusion.GeneNames[FS_UP]);

        boolean threePrimeProm = refFusionData.hasPromiscuousThreeGene(rnaFusion.GeneNames[FS_DOWN]);

        if(fivePrimeProm && threePrimeProm)
        {
            rnaFusion.setKnownType(PROMISCUOUS_BOTH);
        }
        else if(fivePrimeProm)
        {
            rnaFusion.setKnownType(PROMISCUOUS_5.toString());
        }
        else if(threePrimeProm)
        {
            rnaFusion.setKnownType(PROMISCUOUS_3.toString());
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

        if(homologyLength == 0)
        {
            // the insert sequence may explain the difference
            return breakend.getSV().getSvData().insertSequence().length();
        }

        return (homologyLength / 2) + (homologyLength % 2);
    }

    private static final int MAX_BASE_DIFF = 5;

    public static boolean unsplicedPositionMatch(final SvBreakend breakend, int rnaPosition)
    {
        // checks for a position match within the bounds of uncertainty
        int offsetMargin = getHomologyMargin(breakend);

        return abs(breakend.position() - rnaPosition) <= offsetMargin + MAX_BASE_DIFF;
    }

    public boolean isTranscriptBreakendViableForRnaBoundary(
            final BreakendTransData trans, boolean isUpstream, int breakendPosition,
            int rnaPosition, boolean exactRnaPosition)
    {
        // breakend must fall at or before the RNA boundary but not further upstream than the previous splice acceptor

        // if the RNA boundary is at or before the 2nd exon (which has the first splice acceptor), then the breakend can
        // be upstream as far the previous gene or 100K
        final TranscriptData transData = mGeneTransCache.getTranscriptData(trans.gene().geneId(), trans.transName());

        if (transData == null || transData.exons().isEmpty())
            return false;

        int strand = trans.gene().strand();

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
                    if ((rnaPosition == exonData.End)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition > exonData.Start && rnaPosition < nextExonData.Start))
                    {
                        // in which case check whether the breakend is before the next exon's splice acceptor
                        if (nextExonData != null)
                        {
                            return breakendPosition < nextExonData.Start;
                        }

                        // can't take the last exon
                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.Start)
                            || (!exactRnaPosition && prevExonData != null && rnaPosition < exonData.End && rnaPosition > prevExonData.End))
                    {
                        if(prevExonData != null)
                        {
                            return breakendPosition > prevExonData.End;
                        }

                        return false;
                    }
                }
            }
            else
            {
                if((strand == POS_STRAND && rnaPosition <= exonData.Start && exonData.Rank <= 2)
                || (strand == NEG_STRAND && rnaPosition >= exonData.End && exonData.Rank <= 2))
                {
                    int breakendDistance = abs(breakendPosition - rnaPosition);

                    if(breakendDistance > PRE_GENE_PROMOTOR_DISTANCE || trans.hasNegativePrevSpliceAcceptorDistance())
                        return false;
                    else
                        return true;
                }

                if(strand == POS_STRAND)
                {
                    if ((rnaPosition == exonData.Start)
                    || (!exactRnaPosition && prevExonData != null && rnaPosition > prevExonData.Start && rnaPosition < exonData.Start))
                    {
                        if(prevExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition > prevExonData.Start;
                        }

                        return false;
                    }
                }
                else
                {
                    if ((rnaPosition == exonData.End)
                    || (!exactRnaPosition && nextExonData != null && rnaPosition < nextExonData.End && rnaPosition > exonData.End))
                    {
                        if(nextExonData != null)
                        {
                            // after the previous exon's splice acceptor
                            return breakendPosition < nextExonData.Start;
                        }

                        return false;
                    }
                }
            }
        }

        return false;
    }


}
