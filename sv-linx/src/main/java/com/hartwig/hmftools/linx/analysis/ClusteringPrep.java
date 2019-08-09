package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.addSvToChrBreakendMap;
import static com.hartwig.hmftools.linx.cn.LohEvent.CN_DATA_NO_SV;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_OVERLAP;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusteringPrep
{
    public static int MAX_SIMPLE_DUP_DEL_CUTOFF = 5000000;
    public static int MIN_SIMPLE_DUP_DEL_CUTOFF = 100000;
    private static int DEL_DUP_LENGTH_TRIM_COUNT = 5;
    private static int MAX_ARM_COUNT = 41; // excluding the 5 short arms

    private static final Logger LOGGER = LogManager.getLogger(ClusteringPrep.class);

    public static void populateChromosomeBreakendMap(final List<SvVarData> allVariants, final ClusteringState state)
    {
        // add each SV's breakends to a map keyed by chromosome, with the breakends in order of position lowest to highest
        for (final SvVarData var : allVariants)
        {
            addSvToChrBreakendMap(var, state.getChrBreakendMap());
        }

        // set indices
        for(List<SvBreakend> breakendList : state.getChrBreakendMap().values())
        {
            for (int i = 0; i < breakendList.size(); ++i)
            {
                breakendList.get(i).setChrPosIndex(i);
            }
        }
    }

    public static void setSimpleVariantLengths(ClusteringState state)
    {
        long delCutoffLength = 0;
        long dupCutoffLength = 0;

        List<Long> delLengthsList = Lists.newArrayList();
        List<Long> dupLengthsList = Lists.newArrayList();

        int simpleArmCount = 0;

        for (final Map.Entry<String, List<SvBreakend>> entry : state.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();

            // first check for complex events on the arm since these will be skipped
            boolean pArmHasInversions = false;
            boolean qArmHasInversions = false;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.type() != INV)
                    continue;

                if(!pArmHasInversions && (var.arm(true) == CHROMOSOME_ARM_P || var.arm(false) == CHROMOSOME_ARM_P))
                    pArmHasInversions = true;

                if(!qArmHasInversions && (var.arm(true) == CHROMOSOME_ARM_Q || var.arm(false) == CHROMOSOME_ARM_Q))
                    qArmHasInversions = true;

                if(pArmHasInversions && qArmHasInversions)
                    break;
            }

            // skip chromosome altogether
            if(pArmHasInversions && qArmHasInversions)
                continue;

            if(!pArmHasInversions)
                ++simpleArmCount;

            if(!qArmHasInversions)
                ++simpleArmCount;

            for(final SvBreakend breakend : breakendList)
            {
                if(!breakend.usesStart() || !(breakend.getSV().type() == DEL || breakend.getSV().type() == DUP))
                    continue;

                final SvVarData var = breakend.getSV();

                if(pArmHasInversions)
                {
                    if (var.arm(true) == CHROMOSOME_ARM_P || var.arm(false) == CHROMOSOME_ARM_P)
                        continue;
                }

                if(qArmHasInversions)
                {
                    if (var.arm(true) == CHROMOSOME_ARM_Q || var.arm(false) == CHROMOSOME_ARM_Q)
                        continue;
                }

                if(var.type() == DEL)
                    delLengthsList.add(var.length());
                else if(var.type() == DUP)
                    dupLengthsList.add(var.length());
            }

            // LOGGER.debug("sample({}) chr({}) svCount({} delDups({})", sampleId, chromosome, breakendList.size(), armCount);
        }

        int trimCount = (int)round(simpleArmCount / (double)MAX_ARM_COUNT * DEL_DUP_LENGTH_TRIM_COUNT);

        if(delLengthsList.size() > trimCount)
        {
            Collections.sort(delLengthsList);
            int lengthIndex = delLengthsList.size() - trimCount - 1; // 10 items, index 0 - 9, exclude 5 - 9, select 9
            delCutoffLength = delLengthsList.get(lengthIndex);
        }

        if(dupLengthsList.size() > trimCount)
        {
            Collections.sort(dupLengthsList);
            int lengthIndex = dupLengthsList.size() - trimCount - 1;
            dupCutoffLength = dupLengthsList.get(lengthIndex);
        }

        delCutoffLength = min(max(delCutoffLength, MIN_SIMPLE_DUP_DEL_CUTOFF), MAX_SIMPLE_DUP_DEL_CUTOFF);
        dupCutoffLength = min(max(dupCutoffLength, MIN_SIMPLE_DUP_DEL_CUTOFF), MAX_SIMPLE_DUP_DEL_CUTOFF);

        LOGGER.debug("simple dels count({}) cutoff-length({}), dups count({}) cutoff-length({}) simpleArms({}) trimCount({})",
                delLengthsList.size(), delCutoffLength, dupLengthsList.size(), dupCutoffLength, simpleArmCount, trimCount);

        state.setCutoffLengths(delCutoffLength, dupCutoffLength);
    }

    public static void annotateNearestSvData(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // mark each SV's nearest other SV and its relationship - neighbouring or overlapping
        for(Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            int breakendCount = breakendList.size();

            for(int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvVarData var = breakend.getSV();

                final SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                final SvBreakend nextBreakend = (i < breakendCount-1) ? breakendList.get(i + 1) : null;

                long closestDistance = -1;
                if(prevBreakend != null && prevBreakend.getSV() != var)
                {
                    long distance = breakend.position() - prevBreakend.position();
                    closestDistance = distance;
                }

                if(nextBreakend != null && nextBreakend.getSV() != var)
                {
                    long distance = nextBreakend.position() - breakend.position();
                    if(closestDistance < 0 || distance < closestDistance)
                        closestDistance = distance;
                }

                if(closestDistance >= 0 && (var.getNearestSvDistance() == -1 || closestDistance < var.getNearestSvDistance()))
                    var.setNearestSvDistance(closestDistance);

                String relationType = "";
                if((prevBreakend != null && prevBreakend.getSV() == var) || (nextBreakend != null && nextBreakend.getSV() == var))
                    relationType = RELATION_TYPE_NEIGHBOUR;
                else
                    relationType = RELATION_TYPE_OVERLAP;

                var.setNearestSvRelation(relationType);
            }
        }
    }

    protected static void associateBreakendCnEvents(final String sampleId, final ClusteringState state)
    {
        // search for breakends that match LOH and Hom-loss events
        // note that LOH-breakend links are established here and then must be tidied up once the sample is complete

        String currentChromosome = "";
        List<SvBreakend> breakendList = null;

        int missedEvents = 0;

        if(state.getLohEventList() != null && !state.getLohEventList().isEmpty())
        {
            for (final LohEvent lohEvent : state.getLohEventList())
            {
                if (!lohEvent.isSvEvent())
                    continue;

                // use the breakend table to find matching SVs
                if (breakendList == null || !currentChromosome.equals(lohEvent.Chromosome))
                {
                    breakendList = state.getChrBreakendMap().get(lohEvent.Chromosome);
                    currentChromosome = lohEvent.Chromosome;
                }

                if (breakendList == null)
                    continue;

                for (final SvBreakend breakend : breakendList)
                {
                    final SvVarData var = breakend.getSV();

                    if (breakend.orientation() == 1 && var.id() == lohEvent.StartSV)
                    {
                        // check for an INV that the correct end is associated
                        boolean skipInvBreakend = var.type() == INV
                                && abs(breakend.position() - lohEvent.PosStart) > abs(breakend.getOtherBreakend().position() - lohEvent.PosStart);

                        if(!skipInvBreakend)
                        {
                            lohEvent.setBreakend(breakend, true);
                            var.getCluster().addLohEvent(lohEvent);
                        }
                    }

                    if (breakend.orientation() == -1 && var.id() == lohEvent.EndSV)
                    {
                        boolean skipInvBreakend = var.type() == INV
                                && abs(breakend.position() - lohEvent.PosEnd) > abs(breakend.getOtherBreakend().position() - lohEvent.PosEnd);

                        if(!skipInvBreakend)
                        {
                            lohEvent.setBreakend(breakend, false);
                            var.getCluster().addLohEvent(lohEvent);
                        }
                    }

                    if (lohEvent.matchedBothSVs())
                        break;
                }

                if (lohEvent.StartSV != CN_DATA_NO_SV && lohEvent.getBreakend(true) == null)
                    ++missedEvents;

                if (lohEvent.EndSV != CN_DATA_NO_SV && lohEvent.getBreakend(false) == null)
                    ++missedEvents;
            }
        }

        if(state.getHomLossList() != null && !state.getHomLossList().isEmpty())
        {
            for (HomLossEvent homLossEvent : state.getHomLossList())
            {
                if (homLossEvent.StartSV == CN_DATA_NO_SV && homLossEvent.EndSV == CN_DATA_NO_SV)
                    continue;

                breakendList = state.getChrBreakendMap().get(homLossEvent.Chromosome);

                if (breakendList == null)
                    continue;

                for (final SvBreakend breakend : breakendList)
                {
                    if (breakend.orientation() == 1 && breakend.getSV().id() == homLossEvent.StartSV)
                    {
                        homLossEvent.setBreakend(breakend, true);
                    }

                    if (breakend.orientation() == -1 && breakend.getSV().id() == homLossEvent.EndSV)
                    {
                        homLossEvent.setBreakend(breakend, false);
                    }

                    if (homLossEvent.matchedBothSVs())
                        break;
                }

                if (homLossEvent.StartSV != CN_DATA_NO_SV && homLossEvent.getBreakend(true) == null)
                    ++missedEvents;

                if (homLossEvent.EndSV != CN_DATA_NO_SV && homLossEvent.getBreakend(false) == null)
                    ++missedEvents;
            }
        }

        if(missedEvents > 0)
        {
            LOGGER.warn("sample({}) missed {} links to LOH and hom-loss events", sampleId, missedEvents);
        }
    }

}
