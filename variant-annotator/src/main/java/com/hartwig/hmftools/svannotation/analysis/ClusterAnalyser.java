package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.DOUBLE_MINUTE;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.DOUBLE_STRANDED_BREAK;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.REPLICATION_EVENT;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_ENCLOSED;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_NEIGHBOURS;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.SV_GROUP_OVERLAP;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.UNPHASED_EVENTS;

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
        // use a 2-dim matrix to record length for dim 1: SV and din:2 start-end
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

                    if ((i == shortI || j == shortJ) && !(i == shortI && j == shortJ)) {
                        // remove the longer of the lengths where a BE is used twice to measure lengths
                        --tiCount;
                    }
                    else {
                        tiLengthsStr += tiLengthsStr.isEmpty() ? String.valueOf(tiLengths[i][j]) : ";" + tiLengths[i][j];
                    }
                }
            }

            cluster.setTICount(tiCount);
            cluster.setTempInsertLengths(tiLengthsStr);
        }


        boolean validReplication = cluster.getConsistencyCount() == 0 && cluster.getSVs().size() - tiCount <= 1;

        if(validReplication)
            cluster.addAnnotation(REPLICATION_EVENT);

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
            neighbourDistance = mUtils.getShortestProximity(var1, var2);
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



}
