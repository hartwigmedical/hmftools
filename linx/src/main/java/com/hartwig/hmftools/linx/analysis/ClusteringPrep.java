package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.SGL_MAPPING_INF;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.addSvToChrBreakendMap;
import static com.hartwig.hmftools.linx.cn.LohEvent.CN_DATA_NO_SV;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.linx.types.LinxConstants.MAX_SIMPLE_DUP_DEL_CUTOFF;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_SIMPLE_DUP_DEL_CUTOFF;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_OVERLAP;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ClusteringPrep
{
    private static final int DEL_DUP_LENGTH_TRIM_COUNT = 5;
    private static final int MAX_ARM_COUNT = 41; // excluding the 5 short arms

    public static void populateChromosomeBreakendMap(final List<SvVarData> allVariants, final ClusteringState state)
    {
        // add each SV's breakends to a map keyed by chromosome, with the breakends in order of position lowest to highest
        for(final SvVarData var : allVariants)
        {
            addSvToChrBreakendMap(var, state.getChrBreakendMap());
        }

        // set indices
        for(List<SvBreakend> breakendList : state.getChrBreakendMap().values())
        {
            for(int i = 0; i < breakendList.size(); ++i)
            {
                breakendList.get(i).setChrPosIndex(i);
            }
        }
    }

    public static void linkSglMappedInferreds(final List<SvVarData> allVariants)
    {
        final List<SvVarData> sgls = allVariants.stream()
                .filter(x -> x.type() == SGL)
                .filter(x -> !x.getSglMappings().isEmpty())
                .collect(Collectors.toList());

        List<SvVarData> infs = allVariants.stream()
                .filter(x -> x.type() == INF)
                .collect(Collectors.toList());

        for(SvVarData sgl : sgls)
        {
            boolean matched = false;

            for(SglMapping mapping : sgl.getSglMappings())
            {
                for(SvVarData inf : infs)
                {
                    if(!mapping.Chromosome.equals(inf.chromosome(true)))
                        continue;

                    if(mapping.Orientation != inf.orientation(true))
                        continue;

                    if(abs(mapping.Position - inf.position(true)) > 1000)
                        continue;

                    // a link has been found
                    final SvVarData newVar = mergeSglMappedInferred(sgl, inf, mapping, sgl.id()); // ++nextSvId

                    LNX_LOGGER.debug("new SV({}) from sgl({}) with mapping to inf({})", newVar.toString(), sgl.posId(), inf.posId());

                    allVariants.add(newVar);
                    newVar.setLinkedSVs(sgl, inf);

                    newVar.addClusterReason(SGL_MAPPING_INF, inf.id());

                    sgl.setLinkedSVs(newVar, inf);
                    inf.setLinkedSVs(newVar, sgl);
                    matched = true;
                    break;
                }

                if(matched)
                    break;
            }
        }
    }

    public static SvVarData mergeSglMappedInferred(final SvVarData sgl, final SvVarData inf, final SglMapping mapping, final int svId)
    {
        final StructuralVariantData sglSvData = sgl.getSvData();
        final StructuralVariantData infSvData = inf.getSvData();

        StructuralVariantType newType;

        if(sglSvData.startChromosome().equals(infSvData.startChromosome()))
        {
            if(sglSvData.startOrientation() == infSvData.startOrientation())
            {
                newType = INV;
            }
            else
            {
                if((sglSvData.startOrientation() == POS_ORIENT) == (sglSvData.startPosition() < infSvData.startPosition()))
                    newType = DEL;
                else
                    newType = DUP;
            }
        }
        else
        {
            newType = BND;
        }

        StructuralVariantData newSvData = ImmutableStructuralVariantData.builder()
                .id(svId)
                .vcfId(sglSvData.vcfId() != null ? sglSvData.vcfId() : "")
                .startChromosome(sglSvData.startChromosome())
                .endChromosome(infSvData.startChromosome())
                .startPosition(sglSvData.startPosition())
                .endPosition(mapping.Position)
                .startOrientation(sglSvData.startOrientation())
                .endOrientation(mapping.Orientation)
                .startHomologySequence(sglSvData.startHomologySequence())
                .endHomologySequence(infSvData.startHomologySequence())
                .junctionCopyNumber(sglSvData.junctionCopyNumber())
                .startAF(sglSvData.startAF())
                .endAF(infSvData.startAF())
                .adjustedStartAF(sglSvData.adjustedStartAF())
                .adjustedEndAF(infSvData.adjustedStartAF())
                .adjustedStartCopyNumber(sglSvData.adjustedStartCopyNumber())
                .adjustedEndCopyNumber(infSvData.adjustedStartCopyNumber())
                .adjustedStartCopyNumberChange(sglSvData.adjustedStartCopyNumberChange())
                .adjustedEndCopyNumberChange(infSvData.adjustedStartCopyNumberChange())
                .insertSequence(sglSvData.insertSequence())
                .type(newType)
                .filter(sglSvData.filter())
                .imprecise(sglSvData.imprecise())
                .qualityScore(sglSvData.qualityScore())
                .event(sglSvData.event())
                .startTumorVariantFragmentCount(sglSvData.startTumorVariantFragmentCount())
                .startTumorReferenceFragmentCount(sglSvData.startTumorReferenceFragmentCount())
                .startNormalVariantFragmentCount(sglSvData.startNormalVariantFragmentCount())
                .startNormalReferenceFragmentCount(sglSvData.startNormalReferenceFragmentCount())
                .endTumorVariantFragmentCount(infSvData.startTumorVariantFragmentCount())
                .endTumorReferenceFragmentCount(infSvData.startTumorReferenceFragmentCount())
                .endNormalVariantFragmentCount(infSvData.startNormalVariantFragmentCount())
                .endNormalReferenceFragmentCount(infSvData.startNormalReferenceFragmentCount())
                .startIntervalOffsetStart(sglSvData.startIntervalOffsetStart())
                .startIntervalOffsetEnd(sglSvData.startIntervalOffsetEnd())
                .endIntervalOffsetStart(infSvData.startIntervalOffsetStart())
                .endIntervalOffsetEnd(infSvData.startIntervalOffsetEnd())
                .inexactHomologyOffsetStart(sglSvData.inexactHomologyOffsetStart())
                .inexactHomologyOffsetEnd(infSvData.inexactHomologyOffsetStart())
                .startLinkedBy(sglSvData.startLinkedBy())
                .endLinkedBy(infSvData.startLinkedBy())
                .startRefContext(sglSvData.startRefContext())
                .endRefContext(infSvData.startRefContext())
                .recovered(sglSvData.recovered())
                .recoveryMethod((sglSvData.recoveryMethod()))
                .recoveryFilter(sglSvData.recoveryFilter())
                .insertSequenceAlignments(sglSvData.insertSequenceAlignments())
                .insertSequenceRepeatClass(sglSvData.insertSequenceRepeatClass())
                .insertSequenceRepeatType(sglSvData.insertSequenceRepeatType())
                .insertSequenceRepeatOrientation(sglSvData.insertSequenceRepeatOrientation())
                .insertSequenceRepeatCoverage(sglSvData.insertSequenceRepeatCoverage())
                .startAnchoringSupportDistance(sglSvData.startAnchoringSupportDistance())
                .endAnchoringSupportDistance(infSvData.startAnchoringSupportDistance())
                .ponCount(0)
                .build();

        return new SvVarData(newSvData);
    }

    public static void setSimpleVariantLengths(ClusteringState state)
    {
        int delCutoffLength = 0;
        int dupCutoffLength = 0;

        List<Integer> delLengthsList = Lists.newArrayList();
        List<Integer> dupLengthsList = Lists.newArrayList();

        int simpleArmCount = 0;

        for(final Map.Entry<String, List<SvBreakend>> entry : state.getChrBreakendMap().entrySet())
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

                if(!pArmHasInversions && (var.arm(true) == P_ARM || var.arm(false) == P_ARM))
                    pArmHasInversions = true;

                if(!qArmHasInversions && (var.arm(true) == Q_ARM || var.arm(false) == Q_ARM))
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
                if(!breakend.usesStart() || !(breakend.type() == DEL || breakend.type() == DUP))
                    continue;

                final SvVarData var = breakend.getSV();

                if(pArmHasInversions)
                {
                    if(var.arm(true) == P_ARM || var.arm(false) == P_ARM)
                        continue;
                }

                if(qArmHasInversions)
                {
                    if(var.arm(true) == Q_ARM || var.arm(false) == Q_ARM)
                        continue;
                }

                if(var.type() == DEL)
                    delLengthsList.add(var.length());
                else if(var.type() == DUP)
                    dupLengthsList.add(var.length());
            }

            // LNX_LOGGER.debug("sample({}) chr({}) svCount({} delDups({})", sampleId, chromosome, breakendList.size(), armCount);
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

        LNX_LOGGER.debug("simple dels count({}) cutoff-length({}), dups count({}) cutoff-length({}) simpleArms({}) trimCount({})",
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

                // work out closest distance to breakends before and after if not the same SV
                int closestDistance = -1;
                if(prevBreakend != null && prevBreakend.getSV() != var)
                {
                    int distance = breakend.position() - prevBreakend.position();
                    closestDistance = distance;
                }

                if(nextBreakend != null && nextBreakend.getSV() != var)
                {
                    int distance = nextBreakend.position() - breakend.position();
                    if(closestDistance < 0 || distance < closestDistance)
                        closestDistance = distance;
                }

                if(closestDistance < 0)
                    continue;

                if(var.getNearestSvDistance() >= 0 && closestDistance >= var.getNearestSvDistance())
                    continue;

                // a new closest breakend
                var.setNearestSvDistance(closestDistance);

                if(var.type() == BND || var.isSglBreakend() || var.type() == INS)
                {
                    // cannot overlap for these types
                    var.setNearestSvRelation(RELATION_TYPE_NEIGHBOUR);
                }
                else
                {
                    if((nextBreakend != null && nextBreakend.getSV() == var) || (prevBreakend != null && prevBreakend.getSV() == var))
                        var.setNearestSvRelation(RELATION_TYPE_NEIGHBOUR);
                    else
                        var.setNearestSvRelation(RELATION_TYPE_OVERLAP);
                }
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

        if(state.getLohEventList() != null)
        {
            for(final LohEvent lohEvent : state.getLohEventList())
            {
                if(!lohEvent.isSvEvent())
                    continue;

                // use the breakend table to find matching SVs
                if(breakendList == null || !currentChromosome.equals(lohEvent.Chromosome))
                {
                    breakendList = state.getChrBreakendMap().get(lohEvent.Chromosome);
                    currentChromosome = lohEvent.Chromosome;
                }

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    final SvVarData var = breakend.getSV();

                    if(breakend.orientation() == 1 && var.id() == lohEvent.StartSV)
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

                    if(breakend.orientation() == -1 && var.id() == lohEvent.EndSV)
                    {
                        boolean skipInvBreakend = var.type() == INV
                                && abs(breakend.position() - lohEvent.PosEnd) > abs(breakend.getOtherBreakend().position() - lohEvent.PosEnd);

                        if(!skipInvBreakend)
                        {
                            lohEvent.setBreakend(breakend, false);
                            var.getCluster().addLohEvent(lohEvent);
                        }
                    }

                    if(lohEvent.matchedBothSVs())
                        break;
                }

                if(lohEvent.StartSV != CN_DATA_NO_SV && lohEvent.getBreakend(true) == null)
                    ++missedEvents;

                if(lohEvent.EndSV != CN_DATA_NO_SV && lohEvent.getBreakend(false) == null)
                    ++missedEvents;
            }
        }

        if(state.getHomLossList() != null && !state.getHomLossList().isEmpty())
        {
            for(HomLossEvent homLossEvent : state.getHomLossList())
            {
                if(homLossEvent.StartSV == CN_DATA_NO_SV && homLossEvent.EndSV == CN_DATA_NO_SV)
                    continue;

                breakendList = state.getChrBreakendMap().get(homLossEvent.Chromosome);

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    if(breakend.orientation() == 1 && breakend.getSV().id() == homLossEvent.StartSV)
                    {
                        homLossEvent.setBreakend(breakend, true);
                    }

                    if(breakend.orientation() == -1 && breakend.getSV().id() == homLossEvent.EndSV)
                    {
                        homLossEvent.setBreakend(breakend, false);
                    }

                    if(homLossEvent.matchedBothSVs())
                        break;
                }

                if(homLossEvent.StartSV != CN_DATA_NO_SV && homLossEvent.getBreakend(true) == null)
                    ++missedEvents;

                if(homLossEvent.EndSV != CN_DATA_NO_SV && homLossEvent.getBreakend(false) == null)
                    ++missedEvents;
            }
        }

        if(missedEvents > 0)
        {
            // a common cause for these are filtered breakends (eg low VAF or inferred pairs)
            LNX_LOGGER.debug("sample({}) missed {} links to LOH and hom-loss events", sampleId, missedEvents);
        }
    }

}
