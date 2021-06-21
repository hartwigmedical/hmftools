package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.hasLowJcn;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.hasLineRepeatClass;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.hasPolyAorTMotif;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.haveLinkedAssemblies;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_BE;
import static com.hartwig.hmftools.linx.types.ResolvedType.LOW_VAF;
import static com.hartwig.hmftools.linx.types.ResolvedType.PAIR_INF;
import static com.hartwig.hmftools.linx.types.ResolvedType.SGL_MAPPED_INF;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.linx.chaining.ChainJcnLimits;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvFiltering
{
    private Map<SvVarData,ResolvedType> mExcludedSVs; // SV and exclusion reason eg duplicate breakends

    private Map<SvBreakend,SvBreakend> mNonSglDuplicateBreakends;

    private final ClusteringState mState;

    public SvFiltering(final ClusteringState state)
    {
        mState = state;
        mExcludedSVs = Maps.newHashMap();
        mNonSglDuplicateBreakends = Maps.newHashMap();
    }

    private static final double LOW_VAF_THRESHOLD = 0.05;
    private static final int SHORT_INV_DISTANCE = 100;
    private static final int ISOLATED_BND_DISTANCE = 5000;

    private static final int PERMITTED_SGL_DUP_BE_DISTANCE = 1;
    private static final int PERMITTED_DUP_BE_DISTANCE = 35;

    public void applyFilters()
    {
        filterBreakends();
    }

    private void filterBreakends()
    {
        mExcludedSVs.clear();

        // filter out duplicate breakends and low-CN change support INVs and BNDs
        // breakends aren't actually removed until all chromosomes have been processed so that the indices can be preserved for various tests
        Map<String,Set<SvBreakend>> breakendRemovalMap = Maps.newHashMap();

        for(Map.Entry<String, List<SvBreakend>> entry : mState.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();
            List<SvBreakend> breakendList = entry.getValue();

            int breakendCount = breakendList.size();

            Set<SvBreakend> removalList = breakendRemovalMap.get(chromosome);

            if(removalList == null)
            {
                removalList = Sets.newHashSet();
                breakendRemovalMap.put(chromosome, removalList);
            }

            for(int i = 0; i < breakendCount; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                // filter out merged SGL-INF pairs
                if(var.isSglBreakend() && var.getLinkedSVs() != null)
                {
                    mExcludedSVs.put(var, SGL_MAPPED_INF);
                    removalList.add(breakend);
                    continue;
                }

                // first check for SGLs already marked for removal
                if(var.type() == SGL && isSingleDuplicateBreakend(breakendList.get(i).getSV()))
                {
                    mExcludedSVs.put(var, DUP_BE);
                    removalList.add(breakend);
                    continue;
                }

                if((var.type() == BND || var.type() == SGL) && isIsolatedLowVafBnd(var))
                {
                    LNX_LOGGER.trace("SV({}) filtered low VAF isolated BND or SGL", var.id());
                    removalList.add(breakend);

                    if(var.type() == BND)
                        removeRemoteBreakend(breakend.getOtherBreakend(), breakendRemovalMap);

                    mExcludedSVs.put(var, LOW_VAF);
                    continue;
                }

                if(i >= breakendCount - 1)
                    break;

                SvBreakend nextBreakend = breakendList.get(i + 1);
                SvVarData nextVar = nextBreakend.getSV();

                if(var.type() == INV && nextVar == var && isLowVafInversion(breakend, nextBreakend))
                {
                    LNX_LOGGER.trace("SV({}) filtered low VAF / CN change INV", var.id());
                    removalList.add(breakend);
                    removalList.add(nextBreakend);
                    mExcludedSVs.put(var, LOW_VAF);
                    continue;
                }

                if(var.type() == INF && nextVar.type() == INF && breakend.orientation() != nextBreakend.orientation()
                && ChainJcnLimits.jcnMatch(var, nextVar) && !mExcludedSVs.containsKey(var) && !mExcludedSVs.containsKey(nextVar))
                {
                    LNX_LOGGER.trace("SV({} & {}) filtered pair of ploidy-match INFs", var.id(), nextVar.id());
                    removalList.add(breakend);
                    removalList.add(nextBreakend);
                    mExcludedSVs.put(var, PAIR_INF);
                    mExcludedSVs.put(nextVar, PAIR_INF);

                    /* creating a cluster creates downstream issues - simpler to leave unclustered
                    SvCluster newCluster = new SvCluster(mState.getNextClusterId());
                    newCluster.addVariant(var);
                    newCluster.addVariant(nextVar);
                    newCluster.setResolved(true, PAIR_INF);
                    */

                    continue;
                }

                if(nextVar == var)
                    continue;

                if(breakend.orientation() != nextBreakend.orientation())
                    continue;

                int distance = nextBreakend.position() - breakend.position();

                if(distance <= PERMITTED_SGL_DUP_BE_DISTANCE && (var.type() == SGL || nextVar.type() == SGL))
                {
                    LNX_LOGGER.trace("SV({}) filtered proximate duplicate breakend", var.type() == SGL ? var.id() : nextVar.id());

                    if(var.type() == SGL)
                    {
                        mExcludedSVs.put(var, DUP_BE);
                        removalList.add(breakend);
                    }
                    else if(nextVar.type() == SGL)
                    {
                        mExcludedSVs.put(nextVar, DUP_BE);
                        removalList.add(nextBreakend);
                    }
                }
                else if(distance <= PERMITTED_DUP_BE_DISTANCE && var.type() == nextVar.type() && var.isEquivBreakend())
                {
                    // 2 non-SGL SVs may be duplicates, so check their other ends
                    SvBreakend otherBe = breakend.getOtherBreakend();
                    SvBreakend nextOtherBe = nextBreakend.getOtherBreakend();

                    if(otherBe.chromosome().equals(nextOtherBe.chromosome())
                    && abs(otherBe.position() - nextOtherBe.position()) <= PERMITTED_DUP_BE_DISTANCE)
                    {
                        // remove both of the duplicates breakends now

                        // select the one with assembly if only has as them
                        if((var.getTIAssemblies(true).isEmpty() && !nextVar.getTIAssemblies(true).isEmpty())
                        || (var.getTIAssemblies(false).isEmpty() && !nextVar.getTIAssemblies(false).isEmpty()))
                        {
                            LNX_LOGGER.trace("SV({}) filtered eqv-duplicate breakend", var.id());

                            mExcludedSVs.put(var, DUP_BE);
                            removalList.add(breakend);

                            if(breakend.chromosome().equals(otherBe.chromosome()))
                            {
                                removalList.add(otherBe);
                            }
                            else
                            {
                                removeRemoteBreakend(otherBe, breakendRemovalMap);
                            }
                        }
                        else
                        {
                            LNX_LOGGER.trace("SV({}) filtered eqv-duplicate breakend", nextVar.id());

                            mExcludedSVs.put(nextVar, DUP_BE);
                            removalList.add(nextBreakend);

                            if(nextBreakend.chromosome().equals(nextOtherBe.chromosome()))
                            {
                                removalList.add(nextOtherBe);
                            }
                            else
                            {
                                removeRemoteBreakend(nextOtherBe, breakendRemovalMap);
                            }
                        }
                    }
                }
                else if(!var.isSglBreakend() && !nextVar.isSglBreakend())
                {
                    checkCandidateSpanningVariant(breakend, nextBreakend, removalList, breakendRemovalMap);
                    checkCandidateSpanningVariant(nextBreakend, breakend, removalList, breakendRemovalMap);
                }
            }
        }

        // now remove filtered breakends
        for(Map.Entry<String,Set<SvBreakend>> entry : breakendRemovalMap.entrySet())
        {
            final Set<SvBreakend> removalList = entry.getValue();

            if(removalList.isEmpty())
                continue;

            final List<SvBreakend> breakendList = mState.getChrBreakendMap().get(entry.getKey());
            removalList.stream().forEach(x -> breakendList.remove(x));

            // and reset indices after excluding breakends
            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                breakend.setChrPosIndex(i);
            }
        }
    }

    public void clusterExcludedVariants(List<SvCluster> clusters)
    {
        for(Map.Entry<SvVarData,ResolvedType> excludedSv : mExcludedSVs.entrySet())
        {
            SvVarData var = excludedSv.getKey();
            ResolvedType exclusionReason = excludedSv.getValue();

            if(var.getCluster() != null)
            {
                SvCluster newCluster = var.getCluster();
                if(!clusters.contains(newCluster))
                    clusters.add(newCluster);

                continue;
            }

            SvCluster newCluster = new SvCluster(mState.getNextClusterId());
            newCluster.addVariant(var);
            newCluster.setResolved(true, exclusionReason);
            clusters.add(newCluster);
        }
    }

    private boolean isIsolatedLowVafBnd(final SvVarData var)
    {
        if(!hasLowJcn(var))
            return false;

        // ignore possible LINE insertions
        if(hasPolyAorTMotif(var) || hasLineRepeatClass(var))
            return false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && var.isSglBreakend())
                continue;

            final SvBreakend breakend = var.getBreakend(se);
            final List<SvBreakend> breakendList = mState.getChrBreakendMap().get(breakend.chromosome());

            int i = breakend.getChrPosIndex();

            boolean isolatedDown = (i == 0) ? true
                    : breakend.position() - breakendList.get(i - 1).position() >= ISOLATED_BND_DISTANCE;

            boolean isolatedUp = (i >= breakendList.size() - 1) ? true
                    : breakendList.get(i + 1).position() - breakend.position() >= ISOLATED_BND_DISTANCE;

            if(!isolatedDown || !isolatedUp)
                return false;
        }

        return true;
    }

    private boolean isLowVafInversion(final SvBreakend breakend, final SvBreakend nextBreakend)
    {
        final SvVarData var = breakend.getSV();

        if(nextBreakend.getSV() != var || nextBreakend.position() - breakend.position() > SHORT_INV_DISTANCE)
            return false;

        return hasLowJcn(var) || var.calcVaf(true) < LOW_VAF_THRESHOLD;
    }

    private boolean isSingleDuplicateBreakend(final SvVarData var)
    {
        return var.type() == SGL && !var.isInferredSgl() && var.isEquivBreakend();
    }

    private void checkCandidateSpanningVariant(final SvBreakend breakend, final SvBreakend nextBreakend,
            final Set<SvBreakend> removalList, Map<String, Set<SvBreakend>> breakendRemovalMap)
    {
        // one variant has an insert sequence potentially matching the assembled TI of another
        if(breakend.getSV().getSvData().insertSequence().length() >= MIN_TEMPLATED_INSERTION_LENGTH
        && !nextBreakend.getSV().getTIAssemblies(!nextBreakend.usesStart()).isEmpty())
        {
            int posLimitDown = nextBreakend.position() - 1;
            int posLimitUp = nextBreakend.position() + 1;

            if(breakend.position() >= posLimitDown && breakend.position() <= posLimitUp)
            {
                final SvBreakend otherBreakend = breakend.getOtherBreakend();

                final SvBreakend nextBreakend2 = mNonSglDuplicateBreakends.get(otherBreakend);

                if(nextBreakend2 == null)
                {
                    mNonSglDuplicateBreakends.put(breakend, nextBreakend);
                    return;
                }

                if(haveLinkedAssemblies(nextBreakend.getSV(), nextBreakend2.getSV(), !nextBreakend.usesStart(), !nextBreakend2.usesStart()))
                {
                    final SvVarData var = breakend.getSV();

                    LNX_LOGGER.trace("SV({}) has duplicate breakends with assembled breakends({} & {})",
                            var.id(), nextBreakend, nextBreakend2);

                    removalList.add(breakend);

                    if(var.type() == BND)
                        removeRemoteBreakend(otherBreakend, breakendRemovalMap);
                    else
                        removalList.add(otherBreakend);

                    mExcludedSVs.put(var, DUP_BE);
                }
            }
        }
    }

    private void removeRemoteBreakend(final SvBreakend breakend, Map<String,Set<SvBreakend>> breakendRemovalMap)
    {
        Set<SvBreakend> otherList = breakendRemovalMap.get(breakend.chromosome());
        if(otherList == null)
        {
            otherList = Sets.newHashSet();
            breakendRemovalMap.put(breakend.chromosome(), otherList);
        }

        otherList.add(breakend);
    }
}
