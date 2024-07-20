package com.hartwig.hmftools.esvee.assembly.phase;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DEL_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isLocalAssemblyCandidate;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.PhaseGroupBuildWriter;

public class LocalGroupBuilder extends ThreadTask
{
    private final Queue<JunctionGroup> mJunctionGroups;
    private final AssemblyConfig mConfig;
    private final PhaseGroupBuildWriter mWriter;

    private final Set<PhaseGroup> mPhaseGroupsSets;
    private final int mJunctionGroupCount;

    public LocalGroupBuilder(final AssemblyConfig config, final Queue<JunctionGroup> junctionGroups, final PhaseGroupBuildWriter writer)
    {
        super("LocalPhaseGroups");

        mConfig = config;
        mWriter = writer;
        mJunctionGroups = junctionGroups;

        mJunctionGroupCount = junctionGroups.size();
        mPhaseGroupsSets = Sets.newHashSet();
    }

    public Set<PhaseGroup> phaseGroups()
    {
        return mPhaseGroupsSets;
    }

    private static final int LOG_COUNT = 10000;

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mJunctionGroups.size();
                int processedCount = mJunctionGroupCount - remainingCount;

                mPerfCounter.start();

                ++processedCount;

                JunctionGroup junctionGroup = mJunctionGroups.remove();

                if((processedCount % LOG_COUNT) == 0)
                {
                    SV_LOGGER.debug("processed {} junction groups into {} local phase groups", processedCount, mPhaseGroupsSets.size());
                }

                formLocalPhaseGroups(junctionGroup);

                stopCheckLog(junctionGroup.toString(), mConfig.PerfLogTime);
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all phase tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void formLocalPhaseGroups(final JunctionGroup junctionGroup)
    {
        if(junctionGroup.junctionAssemblies().size() <= 1)
            return;

        List<JunctionAssembly> posJunctionAssemblies = junctionGroup.junctionAssemblies()
                .stream().filter(x -> x.isForwardJunction()).collect(Collectors.toList());

        List<JunctionAssembly> negJunctionAssemblies = junctionGroup.junctionAssemblies()
                .stream().filter(x -> !x.isForwardJunction()).collect(Collectors.toList());

        if(posJunctionAssemblies.isEmpty() || negJunctionAssemblies.isEmpty())
            return;

        Collections.sort(posJunctionAssemblies, Comparator.comparingInt(x -> x.junction().Position));
        Collections.sort(negJunctionAssemblies, Comparator.comparingInt(x -> x.junction().Position));

        int negIndex = 0;
        JunctionAssembly startNegAssembly = negJunctionAssemblies.get(negIndex);

        for(JunctionAssembly posAssembly : posJunctionAssemblies)
        {
            while(startNegAssembly.junction().Position < posAssembly.junction().Position - PROXIMATE_DEL_LENGTH)
            {
                ++negIndex;

                if(negIndex >= negJunctionAssemblies.size())
                    return;

                startNegAssembly = negJunctionAssemblies.get(negIndex);
            }

            // now test from this negAssembly to the positional upper limit
            int nextNegIndex = negIndex;
            JunctionAssembly negAssembly = startNegAssembly;

            while(negAssembly.junction().Position < posAssembly.junction().Position + PROXIMATE_DEL_LENGTH)
            {
                if(posAssembly.phaseGroup() == null || negAssembly.phaseGroup() != posAssembly.phaseGroup())
                {
                    // check local phasing criteria
                    if(isLocalAssemblyCandidate(posAssembly, negAssembly) || isLocalFacingLinkCandidate(posAssembly, negAssembly))
                    {
                        PhaseGroupBuilder.linkToPhaseGroups(
                                posAssembly.phaseGroup(), posAssembly, negAssembly, mPhaseGroupsSets, null,
                                mWriter, "Local");
                    }
                }

                ++nextNegIndex;

                if(nextNegIndex >= negJunctionAssemblies.size())
                    break;

                negAssembly = negJunctionAssemblies.get(nextNegIndex);
            }
        }
    }

    private static boolean isLocalFacingLinkCandidate(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(first.refSideSoftClips().isEmpty() && second.refSideSoftClips().isEmpty())
            return false;

        for(int i = 0; i <= 1; ++i)
        {
            JunctionAssembly assembly = (i == 0) ? first : second;
            JunctionAssembly otherAssembly = (i == 0) ? second : first;

            int otherMinAlignedPosition = otherAssembly.minAlignedPosition();
            int otherMaxAlignedPosition = otherAssembly.maxAlignedPosition();

            for(RefSideSoftClip refSideSoftClip : assembly.refSideSoftClips())
            {
                if(!positionWithin(refSideSoftClip.Position, otherMinAlignedPosition, otherMaxAlignedPosition))
                    continue;

                for(String readId : refSideSoftClip.readIds())
                {
                    if(otherAssembly.support().stream().anyMatch(x -> x.id().equals(readId)))
                        return true;

                    if(otherAssembly.candidateSupport().stream().anyMatch(x -> x.id().equals(readId)))
                        return true;
                }
            }
        }

        return false;
    }
}
