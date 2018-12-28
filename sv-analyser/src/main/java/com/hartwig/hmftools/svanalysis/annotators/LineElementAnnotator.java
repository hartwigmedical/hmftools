package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.SHORT_DB_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.areVariantsLinkedByDistance;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LineElementAnnotator {

    public static String KNOWN_LINE_ELEMENT = "Known";
    public static String IDENTIFIED_LINE_ELEMENT = "Ident";
    public static String NO_LINE_ELEMENT = "None";
    public static String SUSPECTED_LINE_ELEMENT = "Suspect";

    private static String CSV_LE_TYPE_IDENTIFIED = "Identified";

    private List<GenomeRegion> mKnownLineElements;
    private List<GenomeRegion> mIdentifiedLineElements;
    private static int PERMITTED_DISTANCE = 5000;

    public static String POLY_A_MOTIF = "AAAAAAAAAAA";
    public static String POLY_T_MOTIF = "TTTTTTTTTTT";

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public LineElementAnnotator()
    {
        mKnownLineElements = Lists.newArrayList();
        mIdentifiedLineElements = Lists.newArrayList();
    }

    public void loadLineElementsFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line;
            while ((line = fileReader.readLine()) != null) {

                if(line.contains("Chromosome"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < 4)
                    continue;

                final GenomeRegion genomeRegion = GenomeRegionFactory.create(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]));

                if(items[3].equals(CSV_LE_TYPE_IDENTIFIED))
                    mIdentifiedLineElements.add(genomeRegion);
                else
                    mKnownLineElements.add(genomeRegion);

//                LOGGER.debug("loaded line element: chr({}) pos({}-{})",
//                        genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
            }

            LOGGER.debug("loaded {} known line elements, {} identified", mKnownLineElements.size(), mIdentifiedLineElements.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    public String isLineElement(final SvVarData svData, final boolean useStart)
    {
        if(mKnownLineElements.isEmpty() && mIdentifiedLineElements.isEmpty())
            return NO_LINE_ELEMENT;

        for(final GenomeRegion genomeRegion : mKnownLineElements)
        {
            if(!genomeRegion.chromosome().equals(svData.chromosome(useStart)))
                continue;

            // test if the SV falls within the LE +/- a buffer
            if(svData.position(useStart) >= genomeRegion.start() - PERMITTED_DISTANCE
            && svData.position(useStart) <= genomeRegion.end() + PERMITTED_DISTANCE)
            {
                LOGGER.debug("var({}) found in known line element({}.>{})",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return KNOWN_LINE_ELEMENT;
            }
        }

        for(final GenomeRegion genomeRegion : mIdentifiedLineElements)
        {
            if(!genomeRegion.chromosome().equals(svData.chromosome(useStart)))
                continue;

            // test if the SV falls within the LE +/- a buffer
            if(svData.position(useStart) >= genomeRegion.start() - PERMITTED_DISTANCE
                    && svData.position(useStart) <= genomeRegion.end() + PERMITTED_DISTANCE)
            {
                LOGGER.debug("var({}) found in identified line element({}.>{})",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return IDENTIFIED_LINE_ELEMENT;
            }
        }

        return NO_LINE_ELEMENT;
    }

    public static boolean hasPolyAorTMotif(final SvVarData var)
    {
        return var.getSvData().insertSequence().contains(POLY_A_MOTIF) || var.getSvData().insertSequence() .contains(POLY_T_MOTIF);
    }


    public static void markLineCluster(final SvCluster cluster, int proximityLength)
    {
        /* Identify a suspected LINE element if:
           - has 2+ BND within 5kb NOT forming a short DB bases on the LINE arm
                AND at least one SV also within 5kb having poly A/T INS sequence
                AND either the 2 BNDs going to different chromosomes or forming a short DB on their non-line (remote) arm
           - OR at least 1 BND with a remote SGL forming a 30 base DB (ie on the remote arm)
                AND EITHER at least one SV also within 5kb OR the remote SGL having poly A/T INS sequence

           Resolve the cluster as type = Line if:
            -  has a suspected line element
            -  every variant in the cluster is part of a KNOWN line element
        */

        // isSpecificCluster(cluster);

        boolean hasSuspected = false;
        boolean hasPolyAorT = false;
        long knownCount = cluster.getSVs().stream().filter(SvVarData::inLineElement).count();

        final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        for (Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                if (var.getLineElement(breakend.usesStart()).equals(SUSPECTED_LINE_ELEMENT))
                    continue;

                SvVarData polyATVar = hasPolyAorTMotif(var) ? var : null;
                List<String> uniqueBndChromosomes = Lists.newArrayList();
                boolean hasRemoteShortBndDB = false;
                List<SvVarData> linkingBnds = Lists.newArrayList();

                boolean isSuspectGroup = false;

                if (var.type() == BND)
                {
                    uniqueBndChromosomes.add(breakend.chromosome());
                    uniqueBndChromosomes.add(var.chromosome(!breakend.usesStart()));

                    linkingBnds.add(var);

                    // test for a remote SGL in a DB
                    final SvLinkedPair dbPair = var.getDBLink(!breakend.usesStart());

                    if (dbPair != null && dbPair.length() <= SHORT_DB_LENGTH && dbPair.getOtherSV(var).type() == SGL
                    && dbPair.getOtherSV(var).getCluster() == cluster)
                    {
                        hasRemoteShortBndDB = true;
                        final SvVarData sgl = dbPair.getOtherSV(var);
                        uniqueBndChromosomes.add(sgl.chromosome(true));
                        linkingBnds.add(sgl);

                        if (hasPolyAorTMotif(sgl))
                        {
                            polyATVar = sgl;
                        }
                    }
                }

                if(hasRemoteShortBndDB && polyATVar != null)
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

                        if (nextBreakend.position() - breakend.position() > proximityLength)
                            break;

                        if (polyATVar == null && hasPolyAorTMotif(nextBreakend.getSV()))
                        {
                            polyATVar = nextBreakend.getSV();
                        }

                        if (nextBreakend.getSV().type() == BND)
                        {
                            final SvLinkedPair dbPair = nextBreakend.getSV().getDBLink(nextBreakend.usesStart());
                            final SvLinkedPair nextDbPair = prevBreakend.getSV().getDBLink(prevBreakend.usesStart());

                            if (dbPair == null || dbPair != nextDbPair || dbPair.length() > SHORT_DB_LENGTH)
                            {
                                // now check the other chromosome for this BND or whether it forms a short DB
                                final String otherChr = nextBreakend.getSV().chromosome(!nextBreakend.usesStart());

                                if(!uniqueBndChromosomes.contains(otherChr))
                                {
                                    linkingBnds.add(nextBreakend.getSV());
                                    uniqueBndChromosomes.add(otherChr);
                                }
                                else
                                {
                                    final SvLinkedPair remoteDbPair = nextBreakend.getSV().getDBLink(!nextBreakend.usesStart());

                                    if(remoteDbPair != null && remoteDbPair.length() <= SHORT_DB_LENGTH
                                    && remoteDbPair.getOtherSV(nextBreakend.getSV()) == prevBreakend.getSV())
                                    {
                                        linkingBnds.add(nextBreakend.getSV());
                                        hasRemoteShortBndDB = true;
                                    }
                                }
                            }
                        }

                        if ((uniqueBndChromosomes.size() >= 3 || hasRemoteShortBndDB) && polyATVar != null)
                        {
                            isSuspectGroup = true;
                            break;
                        }
                    }
                }

                if(polyATVar != null)
                    hasPolyAorT = true;

                if (!isSuspectGroup)
                    continue;

                final String linkingIdsStr = linkingBnds.stream().map(x -> x.id()).collect(Collectors.toList()).toString();

                LOGGER.debug("cluster({}) lineChr({}) uniqueChr({}) linkingSVs({}) hasRemoteShortDB({}) polyAT SV({})",
                        cluster.id(), breakend.chromosome(), uniqueBndChromosomes.toString(), linkingIdsStr, hasRemoteShortBndDB, polyATVar.id());

                hasSuspected = true;

                // otherwise mark every breakend in this proximity as suspect line
                var.setLineElement(SUSPECTED_LINE_ELEMENT, breakend.usesStart());

                for (int j = i + 1; j < breakendList.size(); ++j)
                {
                    final SvBreakend nextBreakend = breakendList.get(j);

                    if (nextBreakend.position() - breakend.position() > proximityLength)
                        break;

                    nextBreakend.getSV().setLineElement(SUSPECTED_LINE_ELEMENT, nextBreakend.usesStart());
                }
            }
        }

        if(cluster.getUniqueSvCount() == knownCount && cluster.getTypeCount(BND) >= 1)
        {
            LOGGER.debug("cluster({}) marked as line with all known({})", cluster.id(), knownCount);
            cluster.markAsLine();
            return;
        }
        else if(hasSuspected)
        {
            if(LOGGER.isDebugEnabled())
            {
                long suspectLine = cluster.getSVs().stream().
                        filter(x -> (x.getLineElement(true).equals(SUSPECTED_LINE_ELEMENT) || x.getLineElement(true).equals(SUSPECTED_LINE_ELEMENT))).count();

                long polyAorT = cluster.getSVs().stream().filter(x -> hasPolyAorTMotif(x)).count();

                LOGGER.debug("cluster({}) marked as line with suspect line elements", cluster.id(), suspectLine, polyAorT);
            }

            cluster.markAsLine();
        }
        else if(hasPolyAorT)
        {
            int sglCount = cluster.getTypeCount(SGL);

            if(sglCount >= 1 && sglCount <= 2 && sglCount == cluster.getCount())
            {
                LOGGER.debug("cluster({}) marked as line with poly A/T SGL", cluster.id());
                cluster.markAsLine();
            }
            else if(cluster.getTypeCount(INS) == 1 && cluster.getCount() == 1)
            {
                cluster.markAsLine();
            }
        }
    }

    public static void markLineCluster_old(final SvCluster cluster, int proximityLength)
    {
        /* Identify a suspected LINE element if:
            - has 2+ BND within 5kb AND at least one SV also within 5kb having poly A/T INS sequence OR
            - has at least 1 BND with a remote SGL forming a 30 base DB or at least one SV poly A/T INS sequence

           New rules:
           - has 2+ BND within 5kb NOT forming a DB of < 30 bases
                AND at least one SV also within 5kb having poly A/T INS sequence
           - OR at least 1 BND with a remote SGL forming a 30 base DB (ie on the remote arm)
                AND EITHER at least one SV also within 5kb OR the remote SGL having poly A/T INS sequence

           Resolve the cluster as type = Line if:
            -  touch a SUSPECTED or KNOWN line element AND contain at least 1 Poly A/T ins sequence OR
            -  every variant in the cluster is part of a KNOWN line element
         */

        List<SvVarData> polyAorTSvs = Lists.newArrayList();
        List<SvVarData> suspectLineSvs = Lists.newArrayList();
        List<SvVarData> svList = cluster.getSVs();

        List<SvVarData> knownLineSvs = svList.stream()
                .filter(SvVarData::inLineElement)
                .collect(Collectors.toList());

        if(cluster.getUniqueSvCount() == knownLineSvs.size())
        {
            LOGGER.debug("cluster({}) marked as line with all known({})", cluster.id(), knownLineSvs.size());
            cluster.markAsLine();
            return;
        }

        for(int i = 0; i < svList.size(); ++i)
        {
            final SvVarData var = svList.get(i);

            if(hasPolyAorTMotif(var))
                polyAorTSvs.add(var);

            if(var.type() != BND && var.type() != SGL)
                continue;

            // check proximity to another possible line element
            for(int j = i + 1; j < svList.size(); ++j)
            {
                final SvVarData otherVar = svList.get(j);

                if(otherVar.type() != BND && otherVar.type() != SGL)
                    continue;

                if(var.type() == SGL && otherVar.type() == SGL)
                    continue;

                if(!hasPolyAorTMotif(var) && !hasPolyAorTMotif(otherVar))
                    continue;

                // check proximity
                boolean areProximate = false;
                boolean v1LinkedOnStart = false;
                boolean v2LinkedOnStart = false;

                for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    if (be1 == SVI_END && var.isNullBreakend())
                        continue;

                    boolean v1Start = isStart(be1);

                    for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        if (be2 == SVI_END && otherVar.isNullBreakend())
                            continue;

                        boolean v2Start = isStart(be2);

                        if (areVariantsLinkedByDistance(var, v1Start, otherVar, v2Start, proximityLength))
                        {
                            areProximate = true;
                            v1LinkedOnStart = v1Start;
                            v2LinkedOnStart = v2Start;
                            break;
                        }
                    }

                    if(areProximate)
                        break;
                }

                if(!areProximate)
                    continue;

                boolean markAsLine = false;
                boolean v1LineIsStart = false;
                boolean v2LineIsStart = false;

                boolean requireShortDB = !(var.type() == BND && otherVar.type() == BND);

                if(var.type() == BND && otherVar.type() == BND)
                {
                    v1LineIsStart = v1LinkedOnStart;
                    v2LineIsStart = v2LinkedOnStart;
                    markAsLine = true;
                }

                // search for a DB to mark to the non-line end
                boolean dbFound = false;
                for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    if (be1 == SVI_END && var.isNullBreakend())
                        continue;

                    boolean v1Start = isStart(be1);

                    for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        if (be2 == SVI_END && otherVar.isNullBreakend())
                            continue;

                        boolean v2Start = isStart(be2);

                        // require a short DB for SGL + BND
                        if(var.getDBLink(v1Start) != null && var.getDBLink(v1Start) == otherVar.getDBLink(v2Start))
                        {
                            dbFound = true;
                            v1LineIsStart = !v1Start;
                            v2LineIsStart = !v2Start;

                            if(requireShortDB && var.getDBLink(v1Start).length() <= SHORT_DB_LENGTH)
                            {
                                markAsLine = true;
                            }

                            break;
                        }
                    }

                    if(markAsLine)
                        break;
                }

                if(markAsLine)
                {
                    if(!var.isLineElement(v1LineIsStart))
                    {
                        LOGGER.debug("var({}) marked as suspect line element on {} dbFound({})",
                                var.posId(), v1LineIsStart ? "start" : "end", dbFound);
                        var.setLineElement(SUSPECTED_LINE_ELEMENT, v1LineIsStart);
                        suspectLineSvs.add(var);
                    }

                    if(!otherVar.isLineElement(v2LineIsStart))
                    {
                        LOGGER.debug("var({}) marked as suspect line element on {} dbFound({})",
                                otherVar.posId(), v2LineIsStart ? "start" : "end", dbFound);
                        otherVar.setLineElement(SUSPECTED_LINE_ELEMENT, v2LineIsStart);
                        suspectLineSvs.add(otherVar);
                    }
                }
            }
        }

        // now check for any other variant proximate to a suspected line element
        for(final SvVarData var : svList)
        {
            if(var.inLineElement() || suspectLineSvs.contains(var))
                continue;

            for(final SvVarData lineVar : suspectLineSvs)
            {
                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    if (be1 == SVI_END && var.isNullBreakend())
                        continue;

                    boolean v1Start = isStart(be1);

                    for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        if (be2 == SVI_END && lineVar.isNullBreakend())
                            continue;

                        boolean v2Start = isStart(be2);

                        if (areVariantsLinkedByDistance(var, v1Start, lineVar, v2Start, proximityLength))
                        {
                            LOGGER.debug("var({}) proximate to suspect line SV({})", var.posId(), lineVar.posId());
                            var.setLineElement(SUSPECTED_LINE_ELEMENT, v1Start);
                            var.setLineElement(SUSPECTED_LINE_ELEMENT, v1Start);
                        }
                    }
                }
            }
        }

        if(!polyAorTSvs.isEmpty() && (!suspectLineSvs.isEmpty() || !knownLineSvs.isEmpty()))
        {
            LOGGER.debug("cluster({}) marked as line with known({}) suspect({}) polyA/T({})",
                    cluster.id(), knownLineSvs.size(), suspectLineSvs.size(), polyAorTSvs.size());
            cluster.markAsLine();
        }
    }

    @Deprecated
    public void setSuspectedLineElements(final Map<String, List<SvBreakend>> chrBreakendMap, int proximityLength)
    {
        // if there are 2 or more BNDs or a BND and SGL within the standard proximity window
        // and have the Poly A or T motif, then classify these as suspected LINE elements
        for(Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            List<Long> positions = Lists.newArrayList();
            List<SvBreakend> potentialLineSVs = Lists.newArrayList();

            boolean isLineGroup = false;

            for (final SvBreakend breakend : breakendList)
            {
                long newPosition = breakend.position();

                // remove earlier breakends which are now too far from the current one
                int index = 0;
                while(index < positions.size())
                {
                    if(newPosition - positions.get(0) <= proximityLength)
                        break;

                    positions.remove(index);
                    potentialLineSVs.remove(index);
                    isLineGroup = false; // require reassessment below
                }

                positions.add(newPosition);
                potentialLineSVs.add(breakend);

                if(positions.size() < 3)
                {
                    isLineGroup = false;
                    continue;
                }

                if(!isLineGroup)
                {
                    boolean hasMultipleRemoteArms = false;
                    int bndCount = 0;
                    boolean hasPolyATMotify = false;
                    String currentOtherChr = "";

                    for (final SvBreakend lineBreakend : potentialLineSVs)
                    {
                        final SvVarData var = lineBreakend.getSV();

                        if (var.type() == BND)
                        {
                            ++bndCount;

                            if(!hasPolyATMotify)
                            {
                                if (hasPolyAorTMotif(var))
                                    hasPolyATMotify = true;
                            }

                            final String otherChr = var.chromosome(!lineBreakend.usesStart());

                            if (currentOtherChr.isEmpty())
                            {
                                currentOtherChr = otherChr;
                            }
                            else if (currentOtherChr.equals(otherChr))
                            {
                                continue;
                            }
                            else
                            {
                                hasMultipleRemoteArms = true;
                            }
                        }

                        if (bndCount >= 2 && hasMultipleRemoteArms && hasPolyATMotify)
                        {
                            isLineGroup = true;
                            break;
                        }
                    }
                }

                if(isLineGroup)
                {
                    for (SvBreakend lineBreakend : potentialLineSVs)
                    {
                        final SvVarData var = lineBreakend.getSV();
                        if(!var.isLineElement(lineBreakend.usesStart()))
                        {
                            LOGGER.debug("var({}) marked as suspect line element", var.posId());
                            var.setLineElement(SUSPECTED_LINE_ELEMENT, lineBreakend.usesStart());
                        }
                    }
                }
            }
        }

    }

}
