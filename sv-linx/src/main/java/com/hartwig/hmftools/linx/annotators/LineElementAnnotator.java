package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.refGenomeChromosome;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.annotators.LineElementType.KNOWN;
import static com.hartwig.hmftools.linx.annotators.LineElementType.SUSPECT;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_DEL_LENGTH;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.SvRegion;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class LineElementAnnotator {

    private final List<SvRegion> mKnownLineElements;
    private PseudoGeneFinder mPseudoGeneFinder;
    private final int mProximityDistance;

    public static final int LINE_ELEMENT_PROXIMITY_DISTANCE = 5000;

    public static final String POLY_A_MOTIF = "AAAAAAAAAAA";
    public static final String POLY_T_MOTIF = "TTTTTTTTTTT";

    private static final int LE_COL_CHR = 0;
    private static final int LE_COL_POS_START = 1;
    private static final int LE_COL_POS_END = 2;

    public LineElementAnnotator(int proximityDistance)
    {
        mProximityDistance = proximityDistance;
        mPseudoGeneFinder = null;
        mKnownLineElements = Lists.newArrayList();
    }

    public void setPseudoGeneFinder(final PseudoGeneFinder pseudoGeneFinder)
    {
        mPseudoGeneFinder = pseudoGeneFinder;
    }

    public void loadLineElementsFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            for(final String line : fileContents)
            {
                if(line.contains("Chromosome"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < LE_COL_POS_END+1)
                    continue;

                final SvRegion lineRegion = new SvRegion(
                        refGenomeChromosome(items[LE_COL_CHR], RG_VERSION),
                        Integer.parseInt(items[LE_COL_POS_START]),
                        Integer.parseInt(items[LE_COL_POS_END]));

                mKnownLineElements.add(lineRegion);
            }

            LNX_LOGGER.info("loaded {} known line elements from file: {}", mKnownLineElements.size(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    public void setKnownLineElements(final SvVarData svData)
    {
        if(mKnownLineElements.isEmpty())
            return;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            for(final SvRegion lineRegion : mKnownLineElements)
            {
                if(!lineRegion.Chromosome.equals(svData.chromosome(isStart(se))))
                    continue;

                // test if the SV falls within the LE +/- a buffer
                if(positionWithin(
                        svData.position(isStart(se)),
                        lineRegion.start() - LINE_ELEMENT_PROXIMITY_DISTANCE,
                        lineRegion.end() + LINE_ELEMENT_PROXIMITY_DISTANCE))
                {
                    LNX_LOGGER.debug("var({}) found in known line element({} -> {})",
                            svData.posId(), lineRegion.Chromosome, lineRegion.start(), lineRegion.end());

                    svData.addLineElement(KNOWN, isStart(se));
                }
            }
        }
    }

    public static boolean hasPolyAorTMotif(final SvVarData var)
    {
        return var.getSvData().insertSequence().contains(POLY_A_MOTIF) || var.getSvData().insertSequence().contains(POLY_T_MOTIF);
    }

    public void markLineCluster(final SvCluster cluster)
    {
        /* Identify a suspected LINE element if:
           - has 2+ BND within 5K NOT forming a short deletion bridge at the suspected LINE source location
                AND at least one SV also within 5kb having poly A/T INS sequence
                AND either the 2 BNDs going to different chromosomes or forming a short DB at their non-line / insertion location
           - OR at least 1 BND with only 1 remote SGL forming a 30 base DB (ie on the remote arm)
                AND EITHER at least one SV also within 5kb OR the remote SGL having poly A/T INS sequence
           - OR there are 2+ breakends within 5kb which both have a polyA insertion sequence with the same orientation

           Resolve the cluster as type = Line if:
            -  has a suspected line element
            -  every variant in the cluster is part of a KNOWN line element
        */

        if(isFilteredResolvedType(cluster.getResolvedType()))
            return;

        boolean hasSuspected = false;
        boolean hasPolyAorT = false;

        // LineClusterState clusterState = new LineClusterState();

        for (Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                // skip if already marked when handled on its other chromosome (ie for a BND)
                if (var.hasLineElement(SUSPECT, breakend.usesStart()))
                    continue;

                final Set<SvVarData> polyAtSVs = Sets.newHashSet();
                final Set<SvBreakend> polyAtBreakends = Sets.newHashSet();

                if(hasPolyAorTMotif(var))
                {
                    polyAtSVs.add(var);
                    polyAtBreakends.add(breakend);
                }

                Set<String> uniqueBndChromosomes = Sets.newHashSet();
                boolean hasRemoteShortBndDB = false;

                boolean isSuspectGroup = false;

                if (var.type() == BND)
                {
                    uniqueBndChromosomes.add(breakend.chromosome());
                    uniqueBndChromosomes.add(var.chromosome(!breakend.usesStart()));

                    // test for a remote SGL in a short DB
                    final SvLinkedPair dbPair = var.getDBLink(!breakend.usesStart());

                    if (isRemoteInsertionDeletionBridge(dbPair) && dbPair.getOtherSV(var).type() == SGL
                    && dbPair.getOtherSV(var).getCluster() == cluster)
                    {
                        hasRemoteShortBndDB = true;
                        final SvVarData sgl = dbPair.getOtherSV(var);

                        if (hasPolyAorTMotif(sgl))
                        {
                            polyAtSVs.add(sgl);
                        }
                    }
                }

                if(hasRemoteShortBndDB && !polyAtSVs.isEmpty())
                {
                    isSuspectGroup = true;
                }
                else
                {
                    // search for proximate BNDs and/or SVs with poly A/T
                    for (int j = i + 1; j < breakendList.size(); ++j)
                    {
                        final SvBreakend prevBreakend = breakendList.get(j - 1);
                        final SvBreakend nextBreakend = breakendList.get(j);

                        if (nextBreakend.position() - breakend.position() > mProximityDistance)
                            break;

                        final SvVarData nextSV = nextBreakend.getSV();

                        if(hasPolyAorTMotif(nextSV))
                        {
                            polyAtSVs.add(nextSV);
                            polyAtBreakends.add(nextBreakend);
                        }

                        if(polyAtBreakends.stream().filter(x -> x.orientation() == 1).count() >= 2
                        || polyAtBreakends.stream().filter(x -> x.orientation() == -1).count() >= 2)
                        {
                            isSuspectGroup = true;
                            break;
                        }

                        if (nextSV.type() == BND)
                        {
                            final SvLinkedPair dbPair = nextBreakend.getDBLink();
                            final SvLinkedPair nextDbPair = prevBreakend.getDBLink();

                            if (dbPair == null || dbPair != nextDbPair || dbPair.length() > MIN_DEL_LENGTH)
                            {
                                // now check the other chromosome for this BND or whether it forms a short DB
                                final String otherChr = nextSV.chromosome(!nextBreakend.usesStart());

                                if(!uniqueBndChromosomes.contains(otherChr))
                                {
                                    uniqueBndChromosomes.add(otherChr);
                                }
                                else
                                {
                                    final SvLinkedPair remoteDbPair = nextSV.getDBLink(!nextBreakend.usesStart());

                                    if(isRemoteInsertionDeletionBridge(remoteDbPair) && remoteDbPair.getOtherSV(nextSV) == prevBreakend.getSV())
                                    {
                                        hasRemoteShortBndDB = true;
                                    }
                                }
                            }
                        }

                        if ((uniqueBndChromosomes.size() >= 3 || hasRemoteShortBndDB) && !polyAtSVs.isEmpty())
                        {
                            isSuspectGroup = true;
                            break;
                        }
                    }
                }

                if(!polyAtSVs.isEmpty())
                    hasPolyAorT = true;

                if (!isSuspectGroup)
                    continue;

                // gather up all breakends within 5K of this suspect line element
                List<SvBreakend> proximateBreakends = Lists.newArrayList(breakend);

                for (int j = i + 1; j < breakendList.size(); ++j)
                {
                    final SvBreakend nextBreakend = breakendList.get(j);

                    if (abs(nextBreakend.position() - breakend.position()) > mProximityDistance)
                        break;

                    proximateBreakends.add(nextBreakend);
                }

                // and in the reverse direction
                for (int j = i - 1; j >= 0; --j)
                {
                    final SvBreakend prevBreakend = breakendList.get(j);

                    if (abs(breakend.position() - prevBreakend.position()) > mProximityDistance)
                        break;

                    proximateBreakends.add(prevBreakend);
                }

                // check for a proximate DEL matching exon positions which would invalidate this as a LINE cluster
                if(mPseudoGeneFinder != null)
                {
                    final SvBreakend pseudoDel = proximateBreakends.stream()
                            .filter(x -> x.getSV().type() == DEL)
                            .filter(x -> mPseudoGeneFinder.variantMatchesPseudogeneExons(x.getSV()))
                            .findFirst().orElse(null);

                    if(pseudoDel != null)
                    {
                        LNX_LOGGER.debug("cluster({}) proximate DEL({}) matches pseudogene exon boundaries",
                                cluster.id(), pseudoDel.getSV().posId());
                        continue;
                    }
                }

                LNX_LOGGER.debug("cluster({}) lineChr({}) uniqueChr({}) remoteShortDB({}) proxPolyATCount({})",
                        cluster.id(), breakend.chromosome(), uniqueBndChromosomes.toString(),hasRemoteShortBndDB, polyAtSVs.size());

                hasSuspected = true;

                // mark every breakend in this proximity as suspect line
                proximateBreakends.forEach(x -> x.getSV().addLineElement(SUSPECT, x.usesStart()));
            }
        }

        markLineCluster(cluster, hasSuspected, hasPolyAorT);
    }

    private boolean isRemoteInsertionDeletionBridge(final SvLinkedPair dbPair)
    {
        if(dbPair == null || dbPair.length() > MIN_DEL_LENGTH)
            return false;

        if(dbPair.firstBreakend().hasLineElement(KNOWN) || dbPair.secondBreakend().hasLineElement(KNOWN))
        {
            return false;
        }

        // no other breakends can be within the proximity cut off
        final List<SvBreakend> breakendList = dbPair.first().getCluster().getChrBreakendMap().get(dbPair.chromosome());

        final SvBreakend lowerBreakend = dbPair.getBreakend(true);
        int startIndex = lowerBreakend.getClusterChrPosIndex();

        if(startIndex > 0)
        {
            if(lowerBreakend.position() - breakendList.get(startIndex - 1).position() < mProximityDistance)
                return false;
        }

        final SvBreakend upperBreakend = dbPair.getBreakend(false);
        int endIndex = upperBreakend.getClusterChrPosIndex();

        if(endIndex < breakendList.size() - 1)
        {
            if(breakendList.get(endIndex + 1).position() - upperBreakend.position() < mProximityDistance)
                return false;
        }

        return true;
    }

    private void markLineCluster(SvCluster cluster, boolean hasSuspected, boolean hasPolyAorT)
    {
        long knownCount = cluster.getSVs().stream().filter(SvVarData::inLineElement).count();

        if(cluster.getSvCount() == knownCount && cluster.getTypeCount(BND) >= 1)
        {
            LNX_LOGGER.debug("cluster({}) marked as line with all known({})", cluster.id(), knownCount);
            cluster.markAsLine();
            return;
        }
        else if(hasSuspected)
        {
            long svInLineCount = cluster.getSVs().stream().filter(x-> x.inLineElement()).count();

            if(LNX_LOGGER.isDebugEnabled())
            {
                long polyAorT = cluster.getSVs().stream().filter(x -> hasPolyAorTMotif(x)).count();

                long suspectLine = cluster.getSVs().stream()
                        .filter(x -> (x.hasLineElement(SUSPECT, true) || x.hasLineElement(SUSPECT, false))).count();

                LNX_LOGGER.debug("cluster({}) anyLine({}) suspect({}) known({}) polyAT({})",
                        cluster.id(), svInLineCount, suspectLine, knownCount, polyAorT);
            }

            if(cluster.getSvCount() <= 10 || svInLineCount * 2 >= cluster.getSvCount())
            {
                cluster.markAsLine();
            }
        }
        else if(hasPolyAorT)
        {
            int sglCount = cluster.getSglBreakendCount();

            if(sglCount >= 1 && sglCount <= 2 && sglCount == cluster.getSvCount())
            {
                LNX_LOGGER.debug("cluster({}) marked as line with poly A/T SGL", cluster.id());
                cluster.markAsLine();
            }
            else if(cluster.getTypeCount(INS) == 1 && cluster.getSvCount() == 1)
            {
                cluster.markAsLine();
            }
        }

    }


}
