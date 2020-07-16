package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.refGenomeChromosome;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.hasLineInsertMotif;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.hasLinePolyAorTMotif;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.hasLineSourceMotif;
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
    private LineClusterState mLineState;

    public static final int LINE_ELEMENT_PROXIMITY_DISTANCE = 5000;

    public static final String POLY_A_MOTIF = "AAAAAAAAAAA";
    public static final String POLY_T_MOTIF = "TTTTTTTTTTT";
    private static final String REPEAT_CLASS_LINE = "LINE/L1";

    private static final int LE_COL_CHR = 0;
    private static final int LE_COL_POS_START = 1;
    private static final int LE_COL_POS_END = 2;

    public LineElementAnnotator(int proximityDistance)
    {
        mProximityDistance = proximityDistance;
        mPseudoGeneFinder = null;
        mKnownLineElements = Lists.newArrayList();
        mLineState = new LineClusterState(proximityDistance);
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
        for(int se = SE_START; se <= SE_END; ++se)
        {
            final SvBreakend breakend = var.getBreakend(se);

            if(breakend == null)
                continue;

            if(hasLineSourceMotif(breakend) || hasLineInsertMotif(breakend))
                return true;
        }

        return false;
    }

    public void markLineCluster(final SvCluster cluster)
    {
        /* Identify a suspected LINE element if:

        1. There are 2+ breakends within 5kb with poly-A/poly-T tails with expected orientations for a source site

        2. There are 2+ BNDs which are not connected at their remote end to a known LINE site (ie within 5kb) with
            - at least one not forming a short DB (< 30 bases) AND
            - at least one breakend within 5kb having a poly-A tail with expected orientation for a source site

        3. There is at least 1 BND with it’s remote breakend proximity clustered with ONLY 1 single breakend AND forming a short DB AND
            - EITHER at least one breakend also within 5kb OR
            - the remote single breakend having a poly-A/poly-T tail with expected orientation for an insertion site
        */

        if(isFilteredResolvedType(cluster.getResolvedType()))
            return;

        boolean hasSuspected = false;

        for (Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            mLineState.clear();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if(breakend.hasLineElement(SUSPECT))
                    continue; // breakends may have already been marked if proximate to an earlier suspect site

                mLineState.addBreakend(breakend);

                if(!mLineState.isSuspected())
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
                            .filter(x -> x.type() == DEL)
                            .filter(x -> mPseudoGeneFinder.variantMatchesPseudogeneExons(x.getSV()))
                            .findFirst().orElse(null);

                    if(pseudoDel != null)
                    {
                        LNX_LOGGER.debug("cluster({}) proximate DEL({}) matches pseudogene exon boundaries",
                                cluster.id(), pseudoDel.getSV().posId());
                        continue;
                    }
                }

                LNX_LOGGER.debug("cluster({}) chromosome({}) marking {} breakends as suspect, state({})",
                        cluster.id(), breakend.chromosome(), proximateBreakends.size(), mLineState.toString());

                // mark every breakend in this local range as suspect line
                for(SvBreakend proxBreakend : proximateBreakends)
                {
                    if(proxBreakend.hasLineElement(SUSPECT))
                        continue;

                    proxBreakend.getSV().addLineElement(SUSPECT, proxBreakend.usesStart());
                    proxBreakend.getSV().addAnnotation(String.format("SLR=%s", mLineState.suspectReason()));
                }

                hasSuspected = true;
            }
        }

        checkIsLineCluster(cluster, hasSuspected, mLineState.hasInsertBreakends());

        mLineState.clear();
    }

    private void checkIsLineCluster(SvCluster cluster, boolean hasSuspected, boolean hasInsertSites)
    {
        /* Resolve the cluster LINE if:

        1. Every non single or inferred breakend variant in the cluster is part of a KNOWN line element

        2. Has a suspected line element AND
        - the cluster has <=10 variants OR at least 50% of the SVs in the cluster have a known or suspected breakend)

        3. The cluster contains only 2 single breakends forming a short DB bridge with one side having an insert poly-A/poly-T tail

        4. The cluster contains exactly 1 single breakend and 1 inferred breakend AND
        - the single breakend has insertSequenceRepeatClass = ‘LINE/L1’
        - (this condition captures LINE insertion sites where we fail to call the polyA side of the insertion)

        5. The cluster contains only 1 single breakend with an insert poly-A/poly-T tail
        */

        boolean allNonSglInKnown = !(cluster.getSVs().stream()
                .filter(x -> !x.isSglBreakend())
                .anyMatch(x -> !(x.hasLineElement(KNOWN, true) || x.hasLineElement(KNOWN, false))));

        if(allNonSglInKnown && cluster.getTypeCount(BND) >= 1)
        {
            long knownCount = cluster.getSVs().stream()
                    .filter(x -> x.hasLineElement(KNOWN, true) || x.hasLineElement(KNOWN, false)).count();

            LNX_LOGGER.debug("cluster({}) marked as line with all known({})", cluster.id(), knownCount);
            cluster.markAsLine();
            return;
        }

        if(hasSuspected)
        {
            long svInLineCount = cluster.getSVs().stream().filter(x-> x.inLineElement()).count();

            if(LNX_LOGGER.isDebugEnabled())
            {
                long polyAorT = cluster.getSVs().stream().filter(x -> hasPolyAorTMotif(x)).count();

                long suspectLine = cluster.getSVs().stream()
                        .filter(x -> (x.hasLineElement(SUSPECT, true) || x.hasLineElement(SUSPECT, false))).count();

                LNX_LOGGER.debug("cluster({}) anyLine({}) suspect({})) polyAT({})",
                        cluster.id(), svInLineCount, suspectLine, polyAorT);
            }

            if(cluster.getSvCount() <= 10)
            {
                cluster.markAsLine();
                return;
            }
            else
            {
                int nonSglCount = cluster.getSvCount() - cluster.getSglBreakendCount();
                int nonSglLineCount = (int)cluster.getSVs().stream().filter(x -> !x.isSglBreakend()).filter(x -> x.inLineElement()).count();

                if(nonSglLineCount * 2 >= nonSglCount)
                {
                    cluster.markAsLine();
                    return;
                }
            }
        }

        // otherwise a 1 or 2 cluster made up of SGLs or INFs only
        if(cluster.getSvCount() <= 2 && cluster.getSglBreakendCount() == cluster.getSvCount())
        {
            /*
            3. The cluster contains only 2 single breakends forming a short DB bridge with one side having an insert poly-A/poly-T tail

            4. The cluster contains exactly 1 single breakend and 1 inferred breakend AND
            - the single breakend has insertSequenceRepeatClass = ‘LINE/L1’
            - (this condition captures LINE insertion sites where we fail to call the polyA side of the insertion)

            5. The cluster contains only 1 single breakend with an insert poly-A/poly-T tail

             */
            if(cluster.getSvCount() == 2 && hasInsertSites)
            {
                final SvVarData var1 = cluster.getSV(0);
                final SvVarData var2 = cluster.getSV(1);

                if(var1.orientation(true) != var2.orientation(true))
                {
                    LNX_LOGGER.debug("cluster({}) marked as line with SGL pair in DB and with insert poly A/T", cluster.id());
                    cluster.markAsLine();
                    return;
                }
            }

            boolean hasSglLineRepeatClass = cluster.getSVs().stream()
                    .filter(x -> x.type() == SGL)
                    .anyMatch(x -> x.getSvData().insertSequenceRepeatClass().equals(REPEAT_CLASS_LINE));

            if(cluster.getTypeCount(SGL) == 1 && (hasSglLineRepeatClass || hasInsertSites))
            {
                LNX_LOGGER.debug("cluster({}) marked as line with SGL with LINE repeat-class({}) insertPolyAT({})",
                        cluster.id(), hasSglLineRepeatClass, hasInsertSites);

                cluster.markAsLine();
                return;
            }
        }
    }

    public void markLineClusterOld(final SvCluster cluster)
    {
        /* Identify a suspected LINE element if:

        1. There are 2+ breakends within 5kb with poly-A/poly-T tails with expected orientations for a source site

        2. There are 2+ BNDs which are not connected at their remote end to a known LINE site (ie within 5kb) with
            - at least one not forming a short DB (< 30 bases) AND
            - at least one breakend within 5kb having a poly-A tail with expected orientation for a source site

        3. There is at least 1 BND with it’s remote breakend proximity clustered with ONLY 1 single breakend AND forming a short DB AND
            - EITHER at least one breakend also within 5kb OR
            - the remote single breakend having a poly-A/poly-T tail with expected orientation for an insertion site
        */

        if(isFilteredResolvedType(cluster.getResolvedType()))
            return;

        boolean hasSuspected = false;
        boolean hasPolyAorT = false;

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
                            .filter(x -> x.type() == DEL)
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

        checkIsLineCluster(cluster, hasSuspected, hasPolyAorT);
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

}
