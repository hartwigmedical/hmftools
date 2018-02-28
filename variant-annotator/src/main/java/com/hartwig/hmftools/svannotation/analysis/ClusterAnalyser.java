package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.DOUBLE_MINUTE;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.DOUBLE_STRANDED_BREAK;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.REPLICATION_EVENT;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_ENCLOSED;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_NEIGHBOURS;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_OVERLAP;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.UNPHASED_EVENTS;

import java.util.List;
import java.util.Vector;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final SvClusteringConfig config, final SvUtilities utils)
    {
        mConfig = config;
        mUtils = utils;

    }

    public void analyserCluster(SvCluster cluster)
    {
        cluster.setConsistencyCount();

        // record an expression of the types of SVs in this cluster
        String clusterTypes = cluster.getClusterTypesAsString();
        cluster.setDesc(clusterTypes);

        cluster.setChromosomalArmCount();

        LOGGER.debug("cluster({}) desc({}) svCount({}) armCount({}) consistent({} count={})",
                cluster.getClusterId(), cluster.getDesc(), cluster.getSVs().size(),
                cluster.getChromosomalArmCount(), cluster.isConsistent(), cluster.getConsistencyCount());

        // for now
        switch(cluster.getSVs().size())
        {
            case 1:
                analyseSingleGroup(cluster);
                break;

            case 2:
                analyseDoubleGroup(cluster);
                break;

            default:
                // not handled for now
                break;
        }

        // what next



    }

    private void analyseSingleGroup(SvCluster cluster)
    {
        SvClusterData var = cluster.getSVs().get(0);

        if(cluster.getConsistencyCount() == 0) {

            LOGGER.debug("cluster({}:{}) is consistent, with type({})", cluster.getClusterId(), cluster.getDesc(), var.type());

            switch (var.type()) {
                case DEL:
                    cluster.addAnnotation(REPLICATION_EVENT);
                    cluster.addAnnotation(DOUBLE_STRANDED_BREAK);
                    break;

                case DUP:
                    cluster.addAnnotation(REPLICATION_EVENT);
                    break;

                case INS:
                    cluster.addAnnotation(REPLICATION_EVENT);
                    break;

                case INV:
                    cluster.addAnnotation(REPLICATION_EVENT);
                    break;

                case BND:
                    cluster.addAnnotation(REPLICATION_EVENT);
                    cluster.addAnnotation(DOUBLE_STRANDED_BREAK);
                    break;
            }
        }
        else
        {
            LOGGER.debug("cluster({}) is inconsistent, with type({})", cluster.getClusterId(), var.type());
        }

        // anything else?
    }

    private void analyseDoubleGroup(SvCluster cluster) {

        SvClusterData var1 = cluster.getSVs().get(0);
        SvClusterData var2 = cluster.getSVs().get(1);

        // count templated insertions
        // can't use the same end twice to count as a DB, so need to find the shortest combo
        int shortestLen = -1;
        int[][] tiLengths = new int[2][2];
        int shortI = 0;
        int shortJ = 0;
        int tiCount = 0;
        String tiLengthsStr = "";

        for(int i = 0; i < 2; ++i)
        {
            boolean v1Start = (i == 0);
            for(int j = 0; j < 2; ++j)
            {
                boolean v2Start = (j == 0);
                tiLengths[i][j] = -1; // initialised

                if(mUtils.areLinkedSection(var1, var2, v1Start, v2Start))
                {
                    int length = mUtils.getProximity(var1, var2, v1Start, v2Start);

                    LOGGER.debug("{}({}) and {}({}) templated insertion, length({}) using {} and {}",
                            var1.type(), var1.posId(), var2.type(), var2.posId(), length,
                            v1Start ? "start" : "end", v2Start ? "start" : "end");

                    tiLengths[i][j] = length;
                    ++tiCount;

                    if(shortestLen == -1 || length < shortestLen)
                    {
                        shortestLen = length;
                        shortI = i;
                        shortJ = j;
                    }
                }
            }
        }

        if(tiCount > 0) {
            for (int i = 0; i < 2; ++i) {

                for (int j = 0; j < 2; ++j) {

                    if(tiLengths[i][j] == -1)
                        continue;

                    if ((i == shortI || j == shortJ) && !(i == shortI && j == shortJ))
                        --tiCount;
                    else
                        tiLengthsStr += tiLengthsStr.isEmpty() ? String.valueOf(tiLengths[i][j]) : ";" + tiLengths[i][j];
                }
            }

            cluster.setTICount(tiCount);
            cluster.setTempInsertLengths(tiLengthsStr);
        }

//        Vector tempInsertLengths = new Vector();
//
//        for(int i = 0; i < 2; ++i)
//        {
//            boolean v1Start = (i == 0);
//            for(int j = 0; j < 2; ++j)
//            {
//                boolean v2Start = (j == 0);
//
//                if(mUtils.areLinkedSection(var1, var2, v1Start, v2Start))
//                {
//                    int length = mUtils.getProximity(var1, var2, v1Start, v2Start);
//
//                    LOGGER.debug("{}({}) and {}({}) templated insertion, length({}) using {} and {}",
//                            var1.type(), var1.posId(), var2.type(), var2.posId(), length,
//                            v1Start ? "start" : "end", v2Start ? "start" : "end");
//
//                    tempInsertLengths.add(new Integer(length));
//                }
//            }
//        }

        boolean validReplication = cluster.getConsistencyCount() == 0 && cluster.getSVs().size() - tiCount <= 1;

        if(validReplication)
            cluster.addAnnotation(REPLICATION_EVENT);

//        if(!tempInsertLengths.isEmpty()) {
//            LOGGER.debug("var({}:{}) and var({}:{}) have {} templated insertions, valid({})",
//                    var1.posId(), var1.type(), var2.posId(), var2.posId(), tempInsertLengths.size(), validReplication);
//        }

        // now count up possible deletion bridges
        int dbCount = 0;
        String dbLengthsStr = "";

        if(mUtils.areTypePair(var1, var2, StructuralVariantType.INV, StructuralVariantType.INV)
        || mUtils.areTypePair(var1, var2, StructuralVariantType.BND, StructuralVariantType.BND)
        || mUtils.areTypePair(var1, var2, StructuralVariantType.DUP, StructuralVariantType.DEL)) {

            // can't use the same end twice to count as a DB, so need to find the shortest combo
            shortestLen = -1;
            int[][] dbLengths = new int[2][2];
            shortI = 0;
            shortJ = 0;

            for (int i = 0; i < 2; ++i) {
                boolean v1Start = (i == 0);
                for (int j = 0; j < 2; ++j) {
                    boolean v2Start = (j == 0);
                    dbLengths[i][j] = -1; // initialised

                    if (mUtils.areSectionBreak(var1, var2, v1Start, v2Start)) {
                        int length = mUtils.getProximity(var1, var2, v1Start, v2Start);

                        LOGGER.debug("{}({}) and {}({}) deletion bridge, length({}) using {} and {}",
                                var1.type(), var1.posId(), var2.type(), var2.posId(), length,
                                v1Start ? "start" : "end", v2Start ? "start" : "end");

                        dbLengths[i][j] = length;
                        ++dbCount;

                        if (shortestLen == -1 || length < shortestLen) {
                            shortestLen = length;
                            shortI = i;
                            shortJ = j;
                        }
                    }
                }
            }

            if(dbCount > 0) {
                for (int i = 0; i < 2; ++i) {

                    for (int j = 0; j < 2; ++j) {

                        if (dbLengths[i][j] == -1)
                            continue;

                        if ((i == shortI || j == shortJ) && !(i == shortI && j == shortJ))
                            --dbCount;
                        else
                            dbLengthsStr += dbLengthsStr.isEmpty() ? String.valueOf(dbLengths[i][j]) : ";" + dbLengths[i][j];
                    }
                }

                cluster.setDSBCount(dbCount);
                cluster.setDSBLengths(dbLengthsStr);
            }
       }

        boolean validDSB = cluster.getSVs().size() == dbCount;

        if(validDSB)
            cluster.addAnnotation(DOUBLE_STRANDED_BREAK);

        if(validReplication || validDSB)
        {
            LOGGER.debug("var({}) and var({}) replication({} count={} lens={}) DSB({} count={} lens={})",
                    var1.posId(), var2.posId(), validReplication, tiCount, tiLengthsStr, validDSB, dbCount, dbLengthsStr);
        }

        // apply some specific tests
        boolean firstEnclosesSecond = mUtils.isWithin(var1, var2);
        boolean secondEnclosesFirst = mUtils.isWithin(var2, var1);
        boolean areEnclosed = firstEnclosesSecond || secondEnclosesFirst;
        boolean areOverlapping = mUtils.isOverlapping(var1, var2);
        boolean areNeighbouring = !firstEnclosesSecond && !secondEnclosesFirst && !areOverlapping;

        // set type from type pairs
        String overlapType = "";

        int neighbourDistance = 0;

        if(areEnclosed)
        {
            overlapType = SV_GROUP_ENCLOSED;
        }
        else if(areOverlapping)
        {
            overlapType = SV_GROUP_OVERLAP;
        }
        else
        {
            overlapType = SV_GROUP_NEIGHBOURS;
            neighbourDistance = mUtils.getNeighbouringProximity(var1, var2);
        }

        LOGGER.debug("cluster({}) desc({}) enclosed({}) overlap({}) neighbours({} dist={})",
                cluster.getClusterId(), cluster.getDesc(), areEnclosed ? (firstEnclosesSecond ? "2nd" : "1st") : "false",
                areOverlapping, areNeighbouring, neighbourDistance);

        if(mUtils.areTypePair(var1, var2, StructuralVariantType.DEL, StructuralVariantType.DEL))
        {
            // cluster.addAnnotation(UNPHASED_EVENTS);

            if(cluster.isConsistent() && (areEnclosed || areOverlapping))
            {
                LOGGER.debug("var({}:{}) and var({}:{}) invalid overlapping DELs)",
                        var1.posId(), var1.type(), var2.posId(), var2.posId());

                cluster.setIsConsistent(false);
            }
            else if(areNeighbouring)
            {
                // FIXME: likely replication event if bridge is short
                // cluster.addAnnotation(REPLICATION_EVENT);
            }
        }
        else if(mUtils.areTypePair(var1, var2, StructuralVariantType.DEL, StructuralVariantType.DUP))
        {
            if(var1.type() == StructuralVariantType.DEL && firstEnclosesSecond
            || var2.type() == StructuralVariantType.DEL && secondEnclosesFirst)
            {
                cluster.addAnnotation(UNPHASED_EVENTS);
                cluster.addAnnotation(DOUBLE_MINUTE);
            }
//            else
//            {
//                // FIXME: if DUP encloses DEL, likely Replication if bridge is short - ??
//
//                cluster.addAnnotation(REPLICATION_EVENT);
//            }
        }


    }

    private void clusterDoubleGroup_Manual(SvCluster cluster)
    {
        SvClusterData var1 = cluster.getSVs().get(0);
        SvClusterData var2 = cluster.getSVs().get(1);

        /* classifications:
            - consistent or not (should have already been established by for now check again)
            - DSB: only a few pairs can satisfy this condition: DEL enclosing DEP, INV-INV and BND-BND
            - replication event - should mostly/all be templated insertions
            - inserted length if relevant
            - originating length if relevant
            - neighbourhood length for SVs on same CRMS
            - possible double-minute
            - phased / 2 events if known

           pairings handled individually for now:
            - DEL-DEL
            - DEL-DUP
            - INV-INV
            - BDN-BND
            - INV + other
            - BND + other
            - INS not handled for now
        */

        // common tests
        boolean firstEnclosesSecond = mUtils.isWithin(var1, var2);
        boolean secondEnclosesFirst = mUtils.isWithin(var2, var1);
        boolean areEnclosed = firstEnclosesSecond || secondEnclosesFirst;
        boolean areOverlapping = mUtils.isOverlapping(var1, var2);
        boolean areNeighbouring = !firstEnclosesSecond && !secondEnclosesFirst && !areOverlapping;

        // set type from type pairs
        String overlapType = "";

        int neighbourDistance = 0;

        if(areEnclosed)
        {
            overlapType = SV_GROUP_ENCLOSED;
        }
        else if(areOverlapping)
        {
            overlapType = SV_GROUP_OVERLAP;
        }
        else
        {
            overlapType = SV_GROUP_NEIGHBOURS;
            neighbourDistance = mUtils.getNeighbouringProximity(var1, var2);
        }

        if(mUtils.areTypePair(var1, var2, StructuralVariantType.DEL, StructuralVariantType.DEL))
        {
            cluster.addAnnotation(UNPHASED_EVENTS);

            if(areEnclosed || areOverlapping)
            {
                cluster.setIsConsistent(false);
            }
            else if(areNeighbouring)
            {
                // FIXME: likely replication event if bridge is short
                // cluster.addAnnotation(REPLICATION_EVENT);
            }
        }
        else if(mUtils.areTypePair(var1, var2, StructuralVariantType.DEL, StructuralVariantType.DUP))
        {
            if(var1.type() == StructuralVariantType.DEL && firstEnclosesSecond
            || var2.type() == StructuralVariantType.DEL && secondEnclosesFirst)
            {
                cluster.addAnnotation(UNPHASED_EVENTS);
                cluster.addAnnotation(DOUBLE_MINUTE);
            }
            else
            {
                // FIXME: if DUP encloses DEL, likely Replication if bridge is short - ??

                cluster.addAnnotation(REPLICATION_EVENT);
            }
        }
        else
        {
            // all other types involve a INS or BND

            if(!cluster.isConsistent())
            {
                // report on why this pair wasn't consistent
            }

            if(mUtils.areTypePair(var1, var2, StructuralVariantType.INV, StructuralVariantType.INV))
            {



            }
            else if(mUtils.areTypePair(var1, var2, StructuralVariantType.BND, StructuralVariantType.BND))
            {


            }
            else
            {
                // other pairs of an INV or BND with any other non-matching type
            }


        }


        // test for replication path




        // test for double-stranded break



    }

    private void clusterMultipleSvCluster(SvCluster cluster) {
        // count templated insertion bridges and deletion bridges
        List<SvClusterData> svList = Lists.newArrayList(cluster.getSVs());

        int tempInsertCount = 0;

        int index1 = 0;
        while (index1 < svList.size()) {
            SvClusterData var1 = svList.get(index1);

            // compare against all others
            int index2 = 0;
            while (index2 < svList.size()) {
                SvClusterData var2 = svList.get(index2);

                if (mUtils.areLinkedSection(var1, var2, true, true)) {
                    tempInsertCount += 1;
                }
                if (mUtils.areLinkedSection(var1, var2, true, false)) {
                    tempInsertCount += 1;
                }
                if (mUtils.areLinkedSection(var1, var2, false, true)) {
                    tempInsertCount += 1;
                }
                if (mUtils.areLinkedSection(var1, var2, false, false)) {
                    tempInsertCount += 1;
                }

                svList.remove(index1); // index will remain the same and so point to the next item
                svList.remove(index2); // index will remain the same and so point to the next item
            }
        }
    }


    private static void annotateInvalidCluster()
    {
//        SV not called	C-T orientaiton wrong AND/OR no replication path	• SV with high ploidy but no copy number change on 1 or both ends
//        FP SV called	C-T orientaiton wrong AND/OR no replication path	• SV with low ploidy and low copy number change at BOTH ends
//        Related events unclustered	C-T orientaiton wrong AND/OR no replication path	• All SVs have high ploidy and consistent copy number change
//        Unrelated events clustered	No replication path


    }

}
